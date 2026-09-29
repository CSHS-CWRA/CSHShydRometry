# Power law fitted by NLS on the natural scale.

#' Fit rating curve using nls on untransformed data
#'
#' @param q A vector of streamflow data.
#' @param h A vector of stage data.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param wts_code Weighting scheme: `"none"`, `"spec"`, or `"prop"`.
#' @param wts Optional vector of weights when `wts_code = "spec"`.
#' @param wts_tol Convergence tolerance for proportional-weight iteration.
#' @param wts_maxiter Maximum iterations for proportional-weight fitting.
#' @param nls_tol Tolerance for nls convergence.
#' @param nls_maxiter Maximum nls iterations.
#' @return An rc_nls object.
#' @examples
#' fit <- rc_nls(q, h, data = thompson)
#' predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
#'
#' # constant coefficient of variation instead of constant variance
#' fit_prop <- rc_nls(q, h, data = thompson, wts_code = "prop")
#' predict(fit_prop, hpred = c(1, 3, 6), conflev = 0.95)
#' @export
rc_nls <- function(
  q,
  h,
  data = NULL,
  wts_code = c("none", "spec", "prop"),
  wts = NULL,
  wts_tol = 1e-6,
  wts_maxiter = 100,
  nls_tol = 1e-6,
  nls_maxiter = 1000
) {
  ## error checks and warnings
  # q and h may name columns of `data`, or be vectors
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  wts_code <- rlang::arg_match(wts_code)

  ## remove missing observations
  qh <- tidyr::drop_na(data.frame(qobs = q, hobs = h))
  q <- qh$qobs
  h <- qh$hobs

  ## generate starting values using lm on log-transformed data
  cstart <- min(h) - 0.1 * (max(h) - min(h))
  start_lm <- stats::lm(log(q) ~ log(h - cstart))
  astart <- exp(start_lm$coefficients[1])
  bstart <- start_lm$coefficients[2]

  # use nls to determine optimal parameters - no weights or specified weights
  if (wts_code != "prop") {
    if (wts_code == "none") {
      wts <- rep(1, length(q))
    } else {
      wts = wts
    }
    checkmate::assert_numeric(wts, len = length(q))
    mod_nls <- stats::nls(
      q ~ a * (h - c)^b,
      weights = wts,
      start = list(a = astart, b = bstart, c = cstart),
      control = list(tol = nls_tol, maxiter = nls_maxiter)
    )
  } else {
    # fit using proportional weights
    yp <- astart * (h - cstart)^bstart
    wts <- 1 / yp^2
    coefs_old <- c(astart, bstart, cstart)
    for (i in 1:wts_maxiter) {
      mod_nls <- stats::nls(
        q ~ a * (h - c)^b,
        weights = wts,
        start = list(a = coefs_old[1], b = coefs_old[2], c = coefs_old[3]),
        control = list(maxiter = nls_maxiter, tol = nls_tol)
      )
      coefs <- as.numeric(stats::coef(mod_nls))
      max_change <- max(abs((coefs - coefs_old) / coefs_old))
      if (max_change < wts_tol) {
        break
      }
      coefs_old <- coefs
      yp <- stats::predict(mod_nls)
      wts <- 1 / yp^2
    }
  }
  mod_sum <- summary(mod_nls)

  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    qh_obs = qh,
    pars = list(
      a = mod_sum$coefficients[1],
      b = mod_sum$coefficients[2],
      c = mod_sum$coefficients[3]
    ),
    rse = mod_sum$sigma,
    wts_code = wts_code,
    weights = wts,
    model = mod_nls
  )
  structure(outlist, class = c("rc_nls", "rating_curve"))
}


#' Predict method for rc_nls objects
#'
#' @inheritParams predict.rc_log_ols
#' @param object An rc_nls object.
#' @export
predict.rc_nls <- function(
  object,
  ...,
  hpred = NULL,
  conflev = NULL,
  predlev = NULL
) {
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  if (predlim && object$wts_code == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1, finite = TRUE)
  hpred_df <- data.frame(h = hpred)
  mod <- object[["model"]]
  wts_code <- object$wts_code
  yvec <- unname(stats::predict(mod, newdata = hpred_df, ...))
  out_df <- data.frame(h = hpred, fit = yvec)
  if (conflim) {
    if (wts_code == "none" || wts_code == "spec") {
      ci_mat <- investr::predFit(
        mod,
        newdata = hpred_df,
        interval = "confidence",
        level = conflev,
        ...
      )
    } else {
      ci_mat <- nlspw_limits(
        mod,
        type = "confidence",
        level = conflev,
        hpred = hpred
      )
    }
    ci_mat <- ci_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(ci_mat) <- paste0("ci_", colnames(ci_mat))
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  # Under "spec" weights the observation variances are supplied rather than
  # estimated, so a new observation's scatter is not identified. Return the
  # columns as NA rather than dropping them, so that the set of columns
  # depends only on what was asked for (Dan, 2026-08-03; Paul agreed).
  if (predlim && wts_code == "spec") {
    out_df$pi_lwr <- NA_real_
    out_df$pi_upr <- NA_real_
  }
  if (predlim && wts_code != "spec") {
    if (wts_code == "none") {
      pi_mat <- investr::predFit(
        mod,
        newdata = hpred_df,
        interval = "prediction",
        level = predlev,
        ...
      )
    } else {
      pi_mat <- nlspw_limits(
        mod,
        type = "prediction",
        level = predlev,
        hpred = hpred
      )
    }
    pi_mat <- pi_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(pi_mat) <- paste0("pi_", colnames(pi_mat))
    out_df <- cbind(out_df, as.data.frame(pi_mat))
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    out_df <- tibble::as_tibble(out_df)
  }
  out_df
}
