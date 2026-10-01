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
#' @param wts How the scatter of the gaugings is modelled: `wts_none()` (or
#'   `"none"`, the default), `wts_prop()` (or `"prop"`), or `wts_spec()`
#'   with the weights. See [wts]. Under `wts_prop()`
#'   the loess curve is refitted in rounds, like the parametric fits.
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
  wts = wts_none()
) {
  # error checks
  # discharge, stage and wts may use columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  # the weighting scheme; specified weights are evaluated in `data`, and
  # kept aligned with the gaugings that remain
  weighting <- resolve_wts(
    wts,
    data,
    keep = stats::complete.cases(discharge, stage),
    fitter = "rc_loess"
  )
  wts_code <- weighting$type
  wts <- weighting$values
  # remove missing observations
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  # fit a loess curve with the given weights (the reweighting helper passes
  # starting values too, which loess does not need)
  surface <- if (extrapolate) "direct" else "interpolate"
  fit_lo <- function(wts, start = NULL) {
    stats::loess(
      discharge ~ stage,
      data = data.frame(discharge = discharge, stage = stage, wts = wts),
      weights = wts,
      degree = degree,
      span = span,
      control = stats::loess.control(surface = surface)
    )
  }
  irls <- NULL
  if (wts_code == "prop") {
    # proportional weights, by iterative reweighting from the unweighted fit
    unweighted <- fit_lo(rep(1, length(discharge)))
    res <- reweight_in_rounds(
      fit_lo,
      yp = as.numeric(stats::predict(unweighted)),
      start = NULL,
      tol = weighting$tol,
      maxiter = weighting$maxiter
    )
    mod_lo <- res$model
    wts <- res$weights
    irls <- res$irls
  } else {
    if (wts_code == "none") {
      wts <- rep(1, length(discharge))
    }
    mod_lo <- fit_lo(wts)
  }
  qh <- tibble::as_tibble(qh)
  outlist <- list(
    gaugings = qh,
    pars = list(),
    settings = list(
      degree = degree,
      span = span,
      extrapolate = extrapolate,
      wts = weighting
    ),
    weights = wts,
    irls = irls,
    wts = weighting,
    enp = mod_lo$enp,
    rse = mod_lo$s,
    model = mod_lo
  )
  structure(outlist, class = c("rc_loess", "rating_curve"))
}


#' Predict method for rc_loess objects
#'
#' @inheritParams predict.rc_power
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
    stage <- stage_grid(object)
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
  out_df <- tibble::as_tibble(out_df)
  out_df
}
