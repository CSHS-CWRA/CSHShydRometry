# Interval methods that resample or simulate.

#' Bootstrap confidence & prediction limits for a two-segment nls rating curve
#'
#' Case-resampling ("pairs") bootstrap: repeatedly resample the gaugings with
#' replacement, refit the two-segment model on each resample, and summarise the
#' resulting spread of fitted rating curves. Returns the same column layout as
#' [predict.rc_nls_2seg()] so the two can be plotted side by side.
#'
#' Why bootstrap? The delta-method limits in [predict.rc_nls_2seg()] (via
#' [investr::predFit()]) linearise the mean function about the fitted
#' parameters, but the two-segment mean is not differentiable in the breakpoint
#' `k` (it is an `ifelse` at `h = k`). That produces an artificial, near-
#' discontinuous widening of the delta-method band at the transition. The
#' bootstrap makes no smoothness assumption, so comparing the two is a direct
#' check on whether the widening near `k` is real or a delta-method artifact.
#'
#' Each resample is refit with the arguments recorded in `object$fit_args`, so
#' the resamples reproduce the original call by construction rather than
#' relying on the caller to restate it.
#'
#' @param object An `rc_nls_2seg` fit (from [rc_nls_2seg()]). The gaugings are
#'   taken from `object$qh_obs` and the fitting arguments from
#'   `object$fit_args`. Under `wts_code = "spec"` the supplied weights are
#'   resampled along with the cases.
#' @param hpred Stage values at which to return limits. Defaults to
#'   [rc_hpred_grid()]: 1000 points spanning the observed stage range.
#' @param ... Currently unused; present so that every `*_limits_2seg()`
#'   function takes the same arguments.
#' @param conflev Coverage for the confidence (mean-curve) interval, or `NULL`.
#' @param predlev Coverage for the prediction (new-observation) interval, or
#'   `NULL`.
#' @param B Number of resamples that must fit successfully.
#' @param seed Optional RNG seed for reproducibility.
#' @param max_tries_factor Cap on total resample attempts (`B * factor`) so a
#'   run terminates even if some resamples fail to converge.
#'
#' @details
#' Confidence limits are the pointwise percentiles of the bootstrap curves
#' (fully nonparametric; this is the part that reveals the delta-method
#' artifact). Prediction limits combine the bootstrap curve uncertainty (the sd
#' of the bootstrap curves at each stage) with the observation-noise sd in
#' quadrature, using a t-quantile on the fit's residual degrees of freedom. The
#' observation-noise sd follows the error model:
#' \itemize{
#'   \item `"none"`: homoscedastic, `sd(q - fitted)`.
#'   \item `"prop"`: constant coefficient of variation,
#'         `fit * sd((q - fitted)/fitted)`.
#'   \item `"spec"`: a new observation's uncertainty is not identified by the
#'         fit, so `pi_lwr`/`pi_upr` are returned as `NA`, as in every other
#'         method.
#' }
#'
#' @return A data frame (tibble if available) with `h`, the point-estimate
#'   curve `fit`, and the requested `ci_lwr`/`ci_upr` and `pi_lwr`/`pi_upr`.
#'   `attr(, "B_success")` records how many resamples converged.
#' @examples
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   fit <- rc_nls_2seg(Q, H, data = sauze, kstart = 1)
#'   boot_limits_2seg(
#'     fit,
#'     hpred = c(1, 2, 4),
#'     conflev = 0.95,
#'     B = 50,
#'     seed = 1
#'   )
#' }
#' @export
boot_limits_2seg <- function(
  object,
  hpred = NULL,
  ...,
  conflev = NULL,
  predlev = NULL,
  B = 1000,
  seed = NULL,
  max_tries_factor = 3
) {
  checkmate::assert_class(object, "rc_nls_2seg")
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1L, finite = TRUE)
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
  fit_args <- object$fit_args
  wts_code <- object$wts_code

  mod0 <- object$model
  qh <- object$qh_obs
  qc <- qh$q
  hc <- qh$h
  n <- length(hc)
  wts_full <- fit_args$wts
  if (wts_code == "spec" && (is.null(wts_full) || length(wts_full) != n)) {
    stop(
      "boot_limits_2seg: spec weights must be supplied and aligned with the ",
      "(NA-dropped) data."
    )
  }
  fit_grid <- as.numeric(
    stats::predict(mod0, newdata = data.frame(h = hpred))
  )

  # Residual pool for prediction limits (model-appropriate scaling).
  mu_obs <- as.numeric(stats::predict(mod0, newdata = data.frame(h = hc)))
  resid_pool <- if (wts_code == "prop") {
    (qc - mu_obs) / mu_obs
  } else {
    (qc - mu_obs)
  }

  # Case-resampling loop. Skip (and retry) resamples that fail to converge.
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(hpred))
  nb <- 0L
  tries <- 0L
  max_tries <- B * max_tries_factor
  while (nb < B && tries < max_tries) {
    tries <- tries + 1L
    s <- sample.int(n, n, replace = TRUE)
    args_b <- fit_args
    if (!is.null(wts_full)) {
      args_b$wts <- wts_full[s]
    }
    fb <- tryCatch(
      do.call(rc_nls_2seg, c(list(q = qc[s], h = hc[s]), args_b)),
      error = function(e) NULL
    )
    if (is.null(fb)) {
      next
    }
    yb <- tryCatch(
      as.numeric(stats::predict(fb$model, newdata = data.frame(h = hpred))),
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
      "boot_limits_2seg: only %d of %d resamples converged.",
      nb,
      B
    ))
    boot_mat <- boot_mat[seq_len(nb), , drop = FALSE]
  }

  out <- data.frame(h = hpred, fit = fit_grid)

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
  if (requireNamespace("tibble", quietly = TRUE)) {
    out <- tibble::as_tibble(out)
  }
  out
}


#' Confidence / prediction limits by simulating the parameter distribution
#'
#' The delta method used by `predict(method = "delta")` fixes the breakpoint at
#' `k-hat` when deciding which segment governs a stage, which makes the band
#' discontinuous at the breakpoint. This function instead draws parameter
#' vectors from their asymptotic distribution `N(theta-hat, Sigma)` and pushes
#' each draw through the model, so every draw carries its own breakpoint and the
#' segment in force varies across draws. Limits are quantiles of the resulting
#' fitted values.
#'
#' This is the same object the delta method approximates -- the sampling
#' distribution of `q-hat(h)` implied by `N(theta-hat, Sigma)` -- but evaluated
#' with no linearisation, in either the parameters or the breakpoint. That
#' matters here: the distribution is skewed with heavy tails away from the
#' gaugings, so quantiles and a symmetric `+/- z * se` disagree substantially.
#'
#' @param object An `rc_nls_2seg` fit (from [rc_nls_2seg()]).
#' @param hpred Stage values at which to return limits. Defaults to
#'   [rc_hpred_grid()]: 1000 points spanning the observed stage range.
#' @param conflev,predlev Coverage for the confidence and prediction intervals,
#'   or `NULL` to omit. Prediction limits add resampled observation noise under
#'   the fit's error model; `NA` for `wts_code = "spec"`.
#' @param M Number of parameter draws.
#' @param seed Optional RNG seed.
#' @param space Parameter scale used to propagate the uncertainty, passed
#'   to [rc_param_space()]. `"original"` uses the fit as-is;
#'   `"log"` refits in gap coordinates, where every draw maps to a
#'   valid curve.
#'
#' @return A data frame (tibble if available) with `h`, `fit` (the fitted curve
#'   at `theta-hat`), and the requested `ci_lwr`/`ci_upr`, `pi_lwr`/`pi_upr`.
#'   `attr(, "frac_invalid")` reports the mean fraction of draws that gave a
#'   non-finite discharge (these are dropped): draws can land outside the region
#'   where the model is defined, which is itself a caution about how well
#'   `N(theta-hat, Sigma)` describes the parameters.
#'
#' @section Caveat:
#' The band inherits `Sigma`, and the breakpoint is a non-regular parameter
#' whose sampling distribution is not normal, so `Var(k-hat)` from the
#' information matrix is unreliable. Where a calibrated band matters, prefer
#' [boot_limits_2seg()], which re-estimates from resampled data instead.
#' @examples
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   fit <- rc_nls_2seg(Q, H, data = sauze, kstart = 1)
#'   sim_limits_2seg(
#'     fit,
#'     hpred = c(1, 2, 4),
#'     conflev = 0.95,
#'     M = 200,
#'     seed = 1
#'   )
#' }
#' @export
sim_limits_2seg <- function(
  object,
  hpred = NULL,
  conflev = NULL,
  predlev = NULL,
  M = 20000,
  seed = NULL,
  space = c("log", "original")
) {
  space <- rlang::arg_match(space)
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("sim_limits_2seg requires the MASS package")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  mod <- object$model
  config <- object$configuration
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1L, finite = TRUE)
  sp <- rc_param_space(object, space)
  TH <- MASS::mvrnorm(M, mu = sp$th, Sigma = sp$Sg)
  colnames(TH) <- names(sp$th)

  fit_curve <- as.numeric(stats::predict(mod, newdata = data.frame(h = hpred)))

  # observation-noise pool for prediction limits
  wts_code <- object$wts_code
  if (!is.null(predlev) && wts_code != "spec") {
    mu_obs <- as.numeric(stats::predict(mod))
    res <- object$qh_obs$q - mu_obs
    if (wts_code == "prop") res <- res / mu_obs
  }

  a <- if (!is.null(conflev)) (1 - conflev) / 2 else NA_real_
  ap <- if (!is.null(predlev)) (1 - predlev) / 2 else NA_real_
  n <- length(hpred)
  ci_l <- ci_u <- pi_l <- pi_u <- rep(NA_real_, n)
  frac_bad <- numeric(n)

  for (i in seq_len(n)) {
    hh <- hpred[i]
    qs <- sp$qmat(hh, TH)
    ok <- is.finite(qs)
    frac_bad[i] <- 1 - mean(ok)
    qs <- qs[ok]
    if (!length(qs)) {
      next
    }
    if (!is.null(conflev)) {
      ci <- stats::quantile(qs, c(a, 1 - a), names = FALSE)
      ci_l[i] <- ci[1]
      ci_u[i] <- ci[2]
    }
    if (!is.null(predlev) && wts_code != "spec") {
      r <- sample(res, length(qs), replace = TRUE)
      ys <- if (wts_code == "prop") qs + qs * r else qs + r
      pj <- stats::quantile(ys, c(ap, 1 - ap), names = FALSE)
      pi_l[i] <- pj[1]
      pi_u[i] <- pj[2]
    }
  }

  out_df <- data.frame(h = hpred, fit = fit_curve)
  if (!is.null(conflev)) {
    out_df$ci_lwr <- ci_l
    out_df$ci_upr <- ci_u
  }
  if (!is.null(predlev)) {
    out_df$pi_lwr <- pi_l
    out_df$pi_upr <- pi_u
  }
  attr(out_df, "frac_invalid") <- mean(frac_bad)
  if (requireNamespace("tibble", quietly = TRUE)) {
    out_df <- tibble::as_tibble(out_df)
  }
  out_df
}
