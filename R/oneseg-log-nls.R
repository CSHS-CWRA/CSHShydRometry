# Power law fitted by NLS on the log-log scale.

#' Fit rating curve using nls on log-transformed data
#'
#' @param q A vector of streamflow data.
#' @param h A vector of stage data.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param tol Tolerance for nls convergence.
#' @return An rc_log_nls object.
#' @examples
#' fit <- rc_log_nls(q, h, data = thompson)
#' predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
#' @export
rc_log_nls <- function(q, h, data = NULL, tol = 1e-6) {
  # q = vector of streamflow data
  # h = vector of stage data
  # q and h may name columns of `data`, or be vectors
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  qh <- tidyr::drop_na(data.frame(qobs = q, hobs = h))
  q <- qh$qobs
  h <- qh$hobs
  # starting estimates
  cstart <- min(h) - 0.1 * (max(h) - min(h))
  lm_mod <- stats::lm(log(q) ~ log(h - cstart))
  b0start <- lm_mod$coef[1]
  b1start <- lm_mod$coef[2]
  # fit model, extract parameters and rse
  mod_nls <- stats::nls(
    log(q) ~ b0 + b1 * log(h - c),
    data = data.frame(q = q, h = h),
    start = list(b0 = b0start, b1 = b1start, c = cstart),
    control = list(tol = tol, maxiter = 1000)
  )
  pars <- stats::coefficients(mod_nls)
  rse <- summary(mod_nls)$sigma
  resids <- summary(mod_nls)$residuals
  a <- exp(pars[1])
  b <- pars[2]
  c <- pars[3]
  # bias corrections
  a_nbc <- a * exp(0.5 * rse^2)
  a_dbc <- (a / length(resids)) * sum(exp(resids))
  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    qh_obs = qh,
    pars = list(a = a, a_nbc = a_nbc, a_dbc = a_dbc, b = b, c = c),
    rse = rse,
    model = mod_nls
  )
  structure(outlist, class = c("rc_log_nls", "rating_curve"))
}


#' Predict method for rc_log_nls objects
#'
#' @inheritParams predict.rc_log_ols
#' @param object An rc_log_nls object.
#' @export
predict.rc_log_nls <- function(
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
  yvec <- unname(exp(stats::predict(mod, newdata = hpred_df, ...)))
  out_df <- data.frame(h = hpred, fit = yvec)
  if (conflim) {
    ci_mat <- exp(investr::predFit(
      mod,
      newdata = hpred_df,
      interval = "confidence",
      level = conflev,
      ...
    ))
    ci_mat <- ci_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(ci_mat) <- paste0("ci_", colnames(ci_mat))
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  if (predlim) {
    pi_mat <- exp(investr::predFit(
      mod,
      newdata = hpred_df,
      interval = "prediction",
      level = predlev,
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
