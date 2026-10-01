# Power law fitted by least squares on the log-log scale.

#' Fit a power-law rating curve on the log-log scale
#'
#' Fits the power law \eqn{Q = a (h - c)^b} by least squares on the log-log
#' scale, `log(discharge) ~ log(a) + b * log(stage - c)`. Equal scatter on that
#' scale means the scatter in discharge is proportional to the flow, and the
#' back-transformed curve estimates the geometric mean discharge at each stage
#' (also the median, when the scatter on the log scale is symmetric).
#'
#' When `c` is estimated, the fit is nonlinear in `c` and is made with
#' [stats::nls()]. When `c` is given, the model is linear on the log-log scale
#' and is fitted with [stats::lm()], and the limits from [predict()] carry no
#' uncertainty in `c`. Fixing `c` at the value estimated by a first fit, as in
#' the examples, treats an estimated `c` as known: the curve is the same, but
#' the limits are narrower than they should be, because the uncertainty in `c`
#' is left out.
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge` and `stage` are
#'   evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param c The stage of zero flow, \eqn{c}: `NULL` (the default) to estimate
#'   it, or a known value, below every gauged stage, to hold it fixed.
#' @param control Settings for [stats::nls()], as from [stats::nls.control()].
#'   Used only when `c` is estimated.
#' @return An `rc_power_log` object; see [rating_curve] for its contents. The
#'   back-transformed coefficient `a` estimates the geometric mean discharge;
#'   `a_corrected` holds two bias-corrected versions for the arithmetic mean:
#'   `nbc`, assuming lognormal errors, and `dbc`, Duan's smearing estimate.
#' @examples
#' fit <- rc_power_log(discharge, stage, data = thompson)
#' coef(fit)
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#'
#' # c known, say from a survey of the control
#' rc_power_log(discharge, stage, data = thompson, c = -1.3)
#'
#' # c estimated, then treated as known: the same curve, with narrower limits
#' # that leave out the uncertainty in c
#' fixed <- rc_power_log(discharge, stage, data = thompson, c = fit$pars$c)
#' predict(fixed, stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_power_log <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  c = NULL,
  control = stats::nls.control(maxiter = 1000, tol = 1e-6)
) {
  # discharge and stage may name columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  checkmate::assert_number(c, null.ok = TRUE, finite = TRUE)
  checkmate::assert_list(control, names = "named")
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  c_input <- c
  if (is.null(c)) {
    # starting estimates from a straight line on the log-log scale
    cstart <- min(stage) - 0.1 * (max(stage) - min(stage))
    lm_mod <- stats::lm(log(discharge) ~ log(stage - cstart))
    mod <- stats::nls(
      log(discharge) ~ b0 + b1 * log(stage - c),
      data = data.frame(discharge = discharge, stage = stage),
      start = list(
        b0 = unname(lm_mod$coef[1]),
        b1 = unname(lm_mod$coef[2]),
        c = cstart
      ),
      control = control
    )
    coefs <- stats::coef(mod)
    log_a <- coefs[["b0"]]
    b <- coefs[["b1"]]
    c <- coefs[["c"]]
  } else {
    if (c >= min(stage)) {
      stop("`c` must be below every gauged stage.", call. = FALSE)
    }
    # with c known the model is linear on the log-log scale
    mod <- stats::lm(
      log(discharge) ~ log(stage - c),
      data = data.frame(discharge = discharge, stage = stage)
    )
    log_a <- unname(stats::coef(mod)[1])
    b <- unname(stats::coef(mod)[2])
  }
  rse <- summary(mod)$sigma
  resids <- stats::residuals(mod)
  a <- exp(log_a)
  # bias corrections, for the arithmetic mean
  a_nbc <- a * exp(0.5 * rse^2)
  a_dbc <- (a / length(resids)) * sum(exp(resids))
  outlist <- list(
    gaugings = tibble::as_tibble(qh),
    pars = list(a = a, b = b, c = c),
    a_corrected = c(nbc = a_nbc, dbc = a_dbc),
    settings = list(c = c_input, control = control),
    rse = rse,
    model = mod
  )
  structure(outlist, class = c("rc_power_log", "rating_curve"))
}


#' Predict method for rc_power_log objects
#'
#' @inheritParams predict.rc_power
#' @param object An `rc_power_log` object.
#' @export
predict.rc_power_log <- function(
  object,
  ...,
  stage = NULL,
  conflev = NULL,
  predlev = NULL
) {
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  if (is.null(stage)) {
    stage <- stage_grid(object)
  }
  checkmate::assert_numeric(stage, min.len = 1, finite = TRUE)
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
  # limits on the log scale, back-transformed; with c estimated, investr
  # accounts for its uncertainty, while with c known an lm does the job
  log_limits <- function(interval, level) {
    lims <- if (inherits(mod, "nls")) {
      investr::predFit(
        mod,
        newdata = stage_df,
        interval = interval,
        level = level,
        ...
      )
    } else {
      stats::predict(
        mod,
        newdata = stage_df,
        interval = interval,
        level = level,
        ...
      )
    }
    exp(lims[, c("lwr", "upr"), drop = FALSE])
  }
  out_df <- data.frame(
    stage = stage,
    fit = unname(exp(stats::predict(mod, newdata = stage_df, ...)))
  )
  if (!is.null(conflev)) {
    ci_mat <- log_limits("confidence", conflev)
    colnames(ci_mat) <- c("ci_lwr", "ci_upr")
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  if (!is.null(predlev)) {
    pi_mat <- log_limits("prediction", predlev)
    colnames(pi_mat) <- c("pi_lwr", "pi_upr")
    out_df <- cbind(out_df, as.data.frame(pi_mat))
  }
  tibble::as_tibble(out_df)
}
