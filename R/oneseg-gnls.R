# Power law fitted by GNLS with a varPower structure.

#' Fit rating curve using gnls
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge` and `stage` are
#'   evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param var_type Variance function for gnls weights (default `nlme::varPower()`).
#' @return An `rc_gnls` object; see [rating_curve] for its contents. The
#'   estimated parameters of the variance function are in `var_pars`.
#' @examples
#' fit <- rc_gnls(discharge, stage, data = thompson)
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_gnls <- function(discharge, stage, ..., data = NULL, var_type = nlme::varPower()) {
  ## error checks
  # discharge and stage may name columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  ## remove missing observations
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  ## generate starting values using lm on log-transformed data
  cstart <- min(stage) - 0.1 * (max(stage) - min(stage))
  start_lm <- stats::lm(log(discharge) ~ log(stage - cstart))
  astart <- as.numeric(exp(start_lm$coefficients[1]))
  bstart <- as.numeric(start_lm$coefficients[2])
  mod_gnls <- nlme::gnls(
    discharge ~ a * (stage - c)^b,
    data = data.frame(discharge = discharge, stage = stage),
    weights = var_type,
    control = nlme::gnlsControl(maxIter = 1e5, minScale = 1e-5),
    start = list(a = astart, b = bstart, c = cstart)
  )
  coefs <- as.numeric(stats::coef(mod_gnls))
  var_pars <- stats::coef(mod_gnls$modelStruct$varStruct, unconstrained = FALSE)
  qh <- tibble::as_tibble(qh)
  outlist <- list(
    gaugings = qh,
    pars = list(a = coefs[1], b = coefs[2], c = coefs[3]),
    var_pars = as.list(var_pars),
    settings = list(var_type = var_type),
    rse = summary(mod_gnls)$sigma,
    model = mod_gnls
  )
  structure(outlist, class = c("rc_gnls", "rating_curve"))
}


#' Predict method for rc_gnls objects
#'
#' @inheritParams predict.rc_log_ols
#' @param object An rc_gnls object.
#' @export
predict.rc_gnls <- function(
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
  yvec <- unname(nlraa::predict_gnls(mod, newdata = stage_df, ...))
  out_df <- data.frame(stage = stage, fit = yvec)
  if (conflim) {
    ci_df <- nlraa::predict_gnls(
      mod,
      newdata = stage_df,
      interval = "confidence",
      level = conflev,
      ...
    )
    ci_mat <- ci_df[, c(3, 4), drop = FALSE]
    colnames(ci_mat) <- c("lwr", "upr")
    colnames(ci_mat) <- paste0("ci_", colnames(ci_mat))
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  if (predlim) {
    pi_df <- nlraa::predict_gnls(
      mod,
      newdata = stage_df,
      interval = "prediction",
      level = predlev,
      ...
    )
    pi_mat <- pi_df[, c(3, 4), drop = FALSE]
    colnames(pi_mat) <- c("lwr", "upr")
    colnames(pi_mat) <- paste0("pi_", colnames(pi_mat))
    out_df <- cbind(out_df, as.data.frame(pi_mat))
  }
  out_df <- tibble::as_tibble(out_df)
  out_df
}
