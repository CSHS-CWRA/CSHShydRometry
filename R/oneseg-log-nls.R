# Power law fitted by NLS on the log-log scale.

#' Fit rating curve using nls on log-transformed data
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge` and `stage` are
#'   evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param control Settings for [stats::nls()], as from [stats::nls.control()].
#' @return An `rc_log_nls` object; see [rating_curve] for its contents. The
#'   back-transformed coefficient `a` estimates the median discharge;
#'   `a_corrected` holds two bias-corrected versions for the mean: `nbc`,
#'   assuming lognormal errors, and `dbc`, Duan's smearing estimate.
#' @examples
#' fit <- rc_log_nls(discharge, stage, data = thompson)
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_log_nls <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  control = stats::nls.control(maxiter = 1000, tol = 1e-6)
) {
  # discharge and stage may name columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  checkmate::assert_list(control, names = "named")
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  # starting estimates
  cstart <- min(stage) - 0.1 * (max(stage) - min(stage))
  lm_mod <- stats::lm(log(discharge) ~ log(stage - cstart))
  b0start <- lm_mod$coef[1]
  b1start <- lm_mod$coef[2]
  # fit model, extract parameters and rse
  mod_nls <- stats::nls(
    log(discharge) ~ b0 + b1 * log(stage - c),
    data = data.frame(discharge = discharge, stage = stage),
    start = list(b0 = b0start, b1 = b1start, c = cstart),
    control = control
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
  qh <- tibble::as_tibble(qh)
  outlist <- list(
    gaugings = qh,
    pars = list(a = unname(a), b = unname(b), c = unname(c)),
    a_corrected = c(nbc = unname(a_nbc), dbc = unname(a_dbc)),
    settings = list(control = control),
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
    ci_mat <- exp(investr::predFit(
      mod,
      newdata = stage_df,
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
      newdata = stage_df,
      interval = "prediction",
      level = predlev,
      ...
    ))
    pi_mat <- pi_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(pi_mat) <- paste0("pi_", colnames(pi_mat))
    out_df <- cbind(out_df, as.data.frame(pi_mat))
  }
  out_df <- tibble::as_tibble(out_df)
  out_df
}
