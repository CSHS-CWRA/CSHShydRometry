# Loess rating curve.

#' Fit rating curve using loess smoother
#'
#' @param q A vector of streamflow data.
#' @param h A vector of stage data.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param degree Degree of local polynomials (1 or 2).
#' @param span Smoothing parameter.
#' @param extrapolate Allow extrapolation beyond observed h range.
#' @param wts_code Weighting scheme: `"none"`, `"spec"`, or `"prop"`.
#' @param wts Optional vector of weights when `wts_code = "spec"`.
#' @return An rc_loess object.
#' @examples
#' fit <- rc_loess(q, h, data = thompson)
#' predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
#' @export
rc_loess <- function(
  q,
  h,
  data = NULL,
  degree = 2,
  span = 0.75,
  extrapolate = TRUE,
  wts_code = c("none", "spec", "prop"),
  wts = NULL
) {
  # error checks
  # q and h may name columns of `data`, or be vectors
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  wts_code <- rlang::arg_match(wts_code)
  if (wts_code == "spec" && !is.numeric(wts)) {
    stop("invalid values of wts")
  }
  # remove missing observations
  qh <- tidyr::drop_na(data.frame(qobs = q, hobs = h))
  q <- qh$qobs
  h <- qh$hobs
  # compute weights
  if (wts_code != "prop") {
    # weights equal or specified
    if (wts_code == "none") {
      wts <- rep(1, length(q))
    }
    if (wts_code == "spec" && !is.numeric(wts)) {
      stop("specified wts are not numeric")
    }
    if (length(wts) != length(q)) warning("length of wts != length of q")
  } else {
    # compute proportional weights - start using loess with no weights
    mod_lo <- stats::loess(q ~ h)
    qp <- stats::predict(mod_lo)
    wts <- 1 / qp^2
  }
  # fit model
  if (extrapolate) {
    mod_lo <- stats::loess(
      q ~ h,
      weights = wts,
      degree = degree,
      span = span,
      control = stats::loess.control(surface = "direct")
    )
  } else {
    mod_lo <- stats::loess(q ~ h, weights = wts, degree = degree, span = span)
  }
  mod_sum <- summary(mod_lo)
  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    qh_obs = qh,
    formula = "q ~ h",
    pars = list(degree = degree, span = span, extrapolate = extrapolate),
    df = mod_sum$df,
    wts_code = wts_code,
    weights = wts,
    rse = mod_sum$sigma,
    model = mod_lo
  )
  structure(outlist, class = c("rc_loess", "rating_curve"))
}


#' Predict method for rc_loess objects
#'
#' @inheritParams predict.rc_log_ols
#' @param object An rc_loess object.
#' @export
predict.rc_loess <- function(
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
  if (predlim) {
    message("Note: prediction limits are not implemented for loess models")
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
    lo_pred <- stats::predict(mod, se = TRUE, newdata = hpred_df, ...)
    tc <- stats::qt(0.5 + 0.5 * conflev, lo_pred$df)
    ci_mat <- cbind(
      lwr = lo_pred$fit - tc * lo_pred$se,
      upr = lo_pred$fit + tc * lo_pred$se
    )
    colnames(ci_mat) <- paste0("ci_", colnames(ci_mat))
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  # prediction limits are not implemented for loess; NA, not absent
  if (predlim) {
    out_df$pi_lwr <- NA_real_
    out_df$pi_upr <- NA_real_
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    out_df <- tibble::as_tibble(out_df)
  }
  out_df
}
