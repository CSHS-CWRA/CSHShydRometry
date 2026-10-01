# Power law fitted by OLS on the log-log scale.

#' Fit rating curve using lm and log-log transform
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge` and `stage` are
#'   evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param c_lwr The lower bound for the c parameter.
#' @param c_upr The upper bound for the c parameter.
#' @return An `rc_log_ols` object; see [rating_curve] for its contents. The
#'   back-transformed coefficient `a` estimates the median discharge;
#'   `a_corrected` holds two bias-corrected versions for the mean: `nbc`,
#'   assuming lognormal errors, and `dbc`, Duan's smearing estimate.
#' @examples
#' fit <- rc_log_ols(discharge, stage, data = thompson)
#' fit
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_log_ols <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  c_lwr = min(stage) - 10 * max(stage),
  c_upr = min(stage) - 0.0001
) {
  # discharge and stage may name columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  # function to minimize
  opt_fun <- function(c, discharge, stage) {
    mod <- stats::lm(log(discharge) ~ log(stage - c))
    summary(mod)$sigma
  }
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  mod_opt <- stats::optim(
    par = c_lwr,
    fn = opt_fun,
    discharge = discharge,
    stage = stage,
    method = "Brent",
    lower = c_lwr,
    upper = c_upr
  )
  c <- mod_opt$par
  mod <- stats::lm(log(discharge) ~ log(stage - c))
  rse <- summary(mod)$sigma
  resids <- summary(mod)$residuals
  a <- exp(mod$coef[1])
  b <- mod$coef[2]
  a_nbc <- a * exp(0.5 * rse^2)
  a_dbc <- (a / length(resids)) * sum(exp(resids))
  qh <- tibble::as_tibble(qh)
  outlist <- list(
    gaugings = qh,
    pars = list(a = unname(a), b = unname(b), c = unname(c)),
    a_corrected = c(nbc = unname(a_nbc), dbc = unname(a_dbc)),
    settings = list(c_lwr = c_lwr, c_upr = c_upr),
    rse = rse,
    model = mod
  )
  structure(outlist, class = c("rc_log_ols", "rating_curve"))
}


#' Predict method for rc_log_ols objects
#'
#' @param object An rc_log_ols object.
#' @param stage Stages at which to predict discharge. If `NULL`, a grid of
#'   1000 equally spaced stages spanning the observed range is used.
#' @param ... Additional arguments passed to the inner `stats::predict()` function.
#' @param conflev The confidence level for the confidence limits; if NULL, no
#'   confidence limits are returned.
#' @param predlev The prediction level for the prediction limits; if NULL, no
#'   prediction limits are returned.
#' @return A tibble with the stages (`stage`), the predicted discharges
#'   (`fit`) and, if requested, the lower and upper confidence limits
#'   (`ci_lwr`, `ci_upr`) and prediction limits (`pi_lwr`, `pi_upr`).
#' @export
predict.rc_log_ols <- function(
  object,
  ...,
  stage = NULL,
  conflev = NULL,
  predlev = NULL
) {
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  if (is.null(stage)) {
    stage <- stage_grid(object)
  }
  checkmate::assert_numeric(stage, min.len = 1, finite = TRUE)
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
  yvec <- unname(exp(stats::predict(mod, newdata = stage_df, ...)))
  out_df <- data.frame(stage = stage, fit = yvec)
  if (conflim) {
    ci_mat <- exp(stats::predict(
      mod,
      interval = "confidence",
      level = conflev,
      newdata = stage_df,
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
      newdata = stage_df,
      ...
    ))
    pi_mat <- pi_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(pi_mat) <- paste0("pi_", colnames(pi_mat))
    out_df <- cbind(out_df, as.data.frame(pi_mat))
  }
  out_df <- tibble::as_tibble(out_df)
  out_df
}
