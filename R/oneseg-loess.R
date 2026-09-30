# Loess rating curve.

#' Fit rating curve using loess smoother
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge`, `stage` and `wts`
#'   are evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param degree Degree of local polynomials (1 or 2).
#' @param span Smoothing parameter.
#' @param extrapolate Allow extrapolation beyond the observed stage range.
#' @param wts_code Weighting scheme: `"none"`, `"spec"`, or `"prop"`.
#' @param wts <[`data-masking`][rlang::args_data_masking]> Weights when
#'   `wts_code = "spec"`, one per gauging: a vector, or an expression
#'   evaluated in `data`, such as `1 / uncertainty_sd^2`.
#' @return An `rc_loess` object; see [rating_curve] for its contents. A loess
#'   curve has no parameters, so `pars` is an empty list.
#' @examples
#' fit <- rc_loess(discharge, stage, data = thompson)
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_loess <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  degree = 2,
  span = 0.75,
  extrapolate = TRUE,
  wts_code = c("none", "spec", "prop"),
  wts = NULL
) {
  # error checks
  # discharge, stage and wts may use columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  wts <- rlang::eval_tidy(rlang::enquo(wts), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  wts_code <- rlang::arg_match(wts_code)
  # remove missing observations
  # keep user-supplied weights aligned with the gaugings that remain
  if (length(wts) == length(discharge)) {
    wts <- wts[stats::complete.cases(discharge, stage)]
  }
  qh <- rc_complete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  wts_input <- wts
  # compute weights
  if (wts_code != "prop") {
    # weights equal or specified
    if (wts_code == "none") {
      wts <- rep(1, length(discharge))
    }
    checkmate::assert_numeric(wts, len = length(discharge), .var.name = "wts")
  } else {
    # compute proportional weights - start using loess with no weights
    mod_lo <- stats::loess(discharge ~ stage)
    qp <- stats::predict(mod_lo)
    wts <- 1 / qp^2
  }
  # fit model
  if (extrapolate) {
    mod_lo <- stats::loess(
      discharge ~ stage,
      weights = wts,
      degree = degree,
      span = span,
      control = stats::loess.control(surface = "direct")
    )
  } else {
    mod_lo <- stats::loess(discharge ~ stage, weights = wts, degree = degree, span = span)
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    gaugings = qh,
    pars = list(),
    settings = list(
      degree = degree,
      span = span,
      extrapolate = extrapolate,
      wts_code = wts_code,
      wts = wts_input
    ),
    weights = wts,
    enp = mod_lo$enp,
    rse = mod_lo$s,
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
  stage = NULL,
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
  if (is.null(stage)) {
    stage <- rc_stage_grid(object)
  }
  checkmate::assert_numeric(stage, min.len = 1, finite = TRUE)
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
  yvec <- unname(stats::predict(mod, newdata = stage_df, ...))
  out_df <- data.frame(stage = stage, fit = yvec)
  if (conflim) {
    lo_pred <- stats::predict(mod, se = TRUE, newdata = stage_df, ...)
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
