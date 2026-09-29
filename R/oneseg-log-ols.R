# Power law fitted by OLS on the log-log scale.

#' Fit rating curve using lm and log-log transform
#'
#' @param q A vector of streamflow data.
#' @param h A vector of stage data.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param c_lwr The lower bound for the c parameter.
#' @param c_upr The upper bound for the c parameter.
#' @return An rc_log_ols object.
#' @examples
#' fit <- rc_log_ols(q, h, data = thompson)
#' fit
#' predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
#' @export
rc_log_ols <- function(
  q,
  h,
  data = NULL,
  c_lwr = min(h) - 10 * max(h),
  c_upr = min(h) - 0.0001
) {
  # q and h may name columns of `data`, or be vectors
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  # function to minimize
  opt_fun <- function(c, q, h) {
    mod <- stats::lm(log(q) ~ log(h - c))
    summary(mod)$sigma
  }
  qh <- tidyr::drop_na(data.frame(qobs = q, hobs = h))
  q <- qh$qobs
  h <- qh$hobs
  mod_opt <- stats::optim(
    par = c_lwr,
    fn = opt_fun,
    q = q,
    h = h,
    method = "Brent",
    lower = c_lwr,
    upper = c_upr
  )
  c <- mod_opt$par
  mod <- stats::lm(log(q) ~ log(h - c))
  rse <- summary(mod)$sigma
  resids <- summary(mod)$residuals
  a <- exp(mod$coef[1])
  b <- mod$coef[2]
  a_nbc <- a * exp(0.5 * rse^2)
  a_dbc <- (a / length(resids)) * sum(exp(resids))
  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    qh_obs = qh,
    pars = list(a = a, a_nbc = a_nbc, a_dbc = a_dbc, b = b, c = c),
    rse = rse,
    model = mod
  )
  structure(outlist, class = c("rc_log_ols", "rating_curve"))
}


#' Predict method for rc_log_ols objects
#'
#' @param object An rc_log_ols object.
#' @param hpred A vector of h values to predict on. If NULL, a grid of 1000
#'   equally spaced values between the minimum and maximum of the input h values
#'   is used.
#' @param ... Additional arguments passed to the inner `stats::predict()` function.
#' @param conflev The confidence level for the confidence limits; if NULL, no
#'   confidence limits are returned.
#' @param predlev The prediction level for the prediction limits; if NULL, no
#'   prediction limits are returned.
#' @return A data frame with the predicted values and confidence/prediction
#'   limits, if requested. The data frame has columns for the predicted values
#'   (fit), the lower confidence limit (ci_lwr), the upper confidence limit
#'   (ci_upr), the lower prediction limit (pi_lwr), and the upper prediction
#'   limit (pi_upr).
#' @export
predict.rc_log_ols <- function(
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
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1, finite = TRUE)
  hpred_df <- data.frame(h = hpred)
  mod <- object[["model"]]
  yvec <- unname(exp(stats::predict(mod, new = hpred_df, ...)))
  out_df <- data.frame(h = hpred, fit = yvec)
  if (conflim) {
    ci_mat <- exp(stats::predict(
      mod,
      interval = "confidence",
      level = conflev,
      new = hpred_df,
      ...
    ))
    ci_mat <- ci_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(ci_mat) <- paste0("ci_", colnames(ci_mat))
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  if (predlim) {
    pi_mat <- exp(stats::predict(
      mod,
      interval = "prediction",
      level = predlev,
      new = hpred_df,
      ...
    ))
    pi_mat <- pi_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(pi_mat) <- paste0("pi_", colnames(pi_mat))
    out_df <- cbind(out_df, as.data.frame(pi_mat))
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    out_df <- tibble::as_tibble(out_df)
  }
  out_df
}
