# Power law fitted by NLS on the natural scale.

#' Fit rating curve using nls on untransformed data
#'
#' @param q A vector of streamflow data.
#' @param h A vector of stage data.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param wts_code Weighting scheme: `"none"`, `"spec"`, or `"prop"`.
#' @param wts Optional vector of weights when `wts_code = "spec"`.
#' @param wts_tol Convergence tolerance under `wts_code = "prop"`: the
#'   reweighting stops once no fitted discharge changes by more than this
#'   fraction from one round to the next.
#' @param wts_maxiter Maximum number of reweighting rounds under
#'   `wts_code = "prop"`. Reaching it gives a warning.
#' @param nls_tol Tolerance for nls convergence.
#' @param nls_maxiter Maximum nls iterations.
#' @return An `rc_nls` object; see [rating_curve] for its contents.
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
  ...,
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
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  wts_code <- rlang::arg_match(wts_code)

  ## remove missing observations
  # keep user-supplied weights aligned with the gaugings that remain
  if (length(wts) == length(q)) {
    wts <- wts[stats::complete.cases(q, h)]
  }
  qh <- rc_complete(q, h)
  q <- qh$q
  h <- qh$h

  ## generate starting values using lm on log-transformed data
  cstart <- min(h) - 0.1 * (max(h) - min(h))
  start_lm <- stats::lm(log(q) ~ log(h - cstart))
  astart <- unname(exp(start_lm$coefficients[1]))
  bstart <- unname(start_lm$coefficients[2])
  wts_input <- wts
  irls <- NULL

  # use nls to determine optimal parameters - no weights or specified weights
  if (wts_code != "prop") {
    if (wts_code == "none") {
      wts <- rep(1, length(q))
    }
    checkmate::assert_numeric(wts, len = length(q), .var.name = "wts")
    mod_nls <- stats::nls(
      q ~ a * (h - c)^b,
      data = data.frame(q = q, h = h),
      weights = wts,
      start = list(a = astart, b = bstart, c = cstart),
      control = list(tol = nls_tol, maxiter = nls_maxiter)
    )
  } else {
    # proportional weights, by iterative reweighting from the log-log curve
    fit_fun <- function(wts, start) {
      stats::nls(
        q ~ a * (h - c)^b,
        data = data.frame(q = q, h = h, wts = wts),
        weights = wts,
        start = start,
        control = list(maxiter = nls_maxiter, tol = nls_tol)
      )
    }
    res <- rc_irls(
      fit_fun,
      yp = astart * (h - cstart)^bstart,
      start = list(a = astart, b = bstart, c = cstart),
      wts_tol = wts_tol,
      wts_maxiter = wts_maxiter
    )
    mod_nls <- res$model
    wts <- res$weights
    irls <- res$irls
  }
  mod_sum <- summary(mod_nls)
  coefs <- stats::coef(mod_nls)

  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    qh_obs = qh,
    pars = list(a = coefs[["a"]], b = coefs[["b"]], c = coefs[["c"]]),
    settings = list(
      wts_code = wts_code,
      wts = wts_input,
      wts_tol = wts_tol,
      wts_maxiter = wts_maxiter,
      nls_tol = nls_tol,
      nls_maxiter = nls_maxiter
    ),
    weights = wts,
    irls = irls,
    rse = mod_sum$sigma,
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
  wts_code <- object$settings$wts_code
  if (predlim && wts_code == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1, finite = TRUE)
  hpred_df <- data.frame(h = hpred)
  mod <- object[["model"]]
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
