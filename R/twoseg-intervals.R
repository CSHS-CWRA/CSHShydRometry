# Bootstrap interval method.

#' Bootstrap limits for a two-segment fit
#'
#' The `method = "boot"` back end of [predict.rc_2seg_powerlaw()], where the
#' method and its arguments are documented. Each resample is refitted with
#' the arguments recorded in `object$settings`, so the resamples reproduce
#' the original call by construction.
#'
#' @param object An `rc_2seg_powerlaw` fit.
#' @param ... Must be empty.
#' @param new_stage,conflev,predlev,B,seed,max_tries_factor As in
#'   [predict.rc_2seg_powerlaw()].
#' @return A tibble, with attribute `"B_success"`.
#' @noRd
boot_limits_2seg <- function(
  object,
  ...,
  new_stage = NULL,
  conflev = NULL,
  predlev = NULL,
  B = 1000,
  seed = NULL,
  max_tries_factor = 3
) {
  checkmate::assert_class(object, "rc_2seg_powerlaw")
  stage <- if (is.null(new_stage)) stage_grid(object) else new_stage
  checkmate::assert_numeric(stage, min.len = 1L, finite = TRUE)
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_count(B, positive = TRUE)
  checkmate::assert_number(seed, null.ok = TRUE)
  rlang::check_dots_empty()

  if (!is.null(seed)) {
    set.seed(seed)
  }
  # The arguments the fit was made with. Refitting with anything else would
  # bootstrap a different model from the one being reported on.
  fit_args <- object$settings
  wts_code <- object$settings$variance$type

  mod0 <- object$model
  qh <- object$gaugings
  qc <- qh$discharge
  hc <- qh$stage
  n <- length(hc)
  wts_full <- if (wts_code == "spec") fit_args$variance$values
  if (wts_code == "spec" && (is.null(wts_full) || length(wts_full) != n)) {
    stop(
      "The bootstrap needs the variances given to `var_spec()`, one for each ",
      "gauging the fit kept.",
      call. = FALSE
    )
  }
  fit_grid <- as.numeric(
    stats::predict(mod0, newdata = data.frame(stage = stage))
  )

  # Residual pool for prediction limits (model-appropriate scaling).
  mu_obs <- as.numeric(stats::predict(mod0, newdata = data.frame(stage = hc)))
  resid_pool <- if (wts_code == "prop") {
    (qc - mu_obs) / mu_obs
  } else {
    (qc - mu_obs)
  }

  # Case-resampling loop. Skip (and retry) resamples that fail to converge.
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(stage))
  nb <- 0L
  tries <- 0L
  max_tries <- B * max_tries_factor
  while (nb < B && tries < max_tries) {
    tries <- tries + 1L
    s <- sample.int(n, n, replace = TRUE)
    args_b <- fit_args
    if (wts_code == "spec") {
      args_b$variance <- new_variance("spec", values = wts_full[s])
    }
    # a resample whose reweighting does not converge counts as a failure,
    # like one whose fit errors, rather than warning once per resample
    fb <- tryCatch(
      suppressWarnings(
        do.call(rc_2seg_powerlaw, c(list(discharge = qc[s], stage = hc[s]), args_b))
      ),
      error = function(e) NULL
    )
    if (is.null(fb) || isFALSE(fb$irls$converged)) {
      next
    }
    yb <- tryCatch(
      as.numeric(stats::predict(fb$model, newdata = data.frame(stage = stage))),
      error = function(e) NULL
    )
    if (is.null(yb) || any(!is.finite(yb))) {
      next
    }
    nb <- nb + 1L
    boot_mat[nb, ] <- yb
  }
  if (nb < B) {
    warning(sprintf(
      "Bootstrap: only %d of %d resamples converged.",
      nb,
      B
    ))
    boot_mat <- boot_mat[seq_len(nb), , drop = FALSE]
  }

  out <- data.frame(stage = stage, fit = fit_grid)

  # Confidence limits: percentiles of the bootstrap mean curves.
  if (!is.null(conflev)) {
    a <- (1 - conflev) / 2
    ci <- t(apply(
      boot_mat,
      2,
      stats::quantile,
      probs = c(a, 1 - a),
      names = FALSE,
      na.rm = TRUE
    ))
    out$ci_lwr <- ci[, 1]
    out$ci_upr <- ci[, 2]
  }

  # Prediction limits: combine the bootstrap curve uncertainty (sd of the
  # bootstrap curves at each stage) with the observation-noise sd in quadrature,
  # using a t-quantile on the fit's residual df. (A fully nonparametric version
  # -- resampling a residual per grid point -- is very jagged because the
  # residual pool is small, so we use the smooth variance-combine instead.)
  if (!is.null(predlev)) {
    if (wts_code == "spec") {
      out$pi_lwr <- NA_real_
      out$pi_upr <- NA_real_
    } else {
      tq <- stats::qt(1 - (1 - predlev) / 2, df = stats::df.residual(mod0))
      # curve uncertainty
      se_boot <- apply(boot_mat, 2, stats::sd, na.rm = TRUE)
      s_obs <- if (wts_code == "prop") {
        fit_grid * stats::sd(resid_pool) # constant-CV noise
      } else {
        rep(stats::sd(resid_pool), length(fit_grid)) # homoscedastic
      }
      half <- tq * sqrt(se_boot^2 + s_obs^2)
      out$pi_lwr <- fit_grid - half
      out$pi_upr <- fit_grid + half
    }
  }

  attr(out, "B_success") <- nb
  out <- tibble::as_tibble(out)
  out
}
