# The structure shared by every fit, and the methods that rely on it.

#' Rating curve fits
#'
#' Every `rc_*()` constructor returns a list with class
#' `c("rc_<method>", "rating_curve")`. These elements are part of the
#' interface:
#'
#' \describe{
#'   \item{`gaugings`}{A tibble of the gaugings used, with columns
#'     `discharge` and `stage`, after dropping any with a missing value.}
#'   \item{`curve_parameters`}{The estimated parameters of the curve, as a named list with
#'     one element per parameter type. For a multi-segment curve an element
#'     holds one value per segment (for `k`, per breakpoint); elements may
#'     differ in length where a parameter is fixed by how the segments join. Empty
#'     for [rc_loess()], which has no parameters. [coef()] returns the same
#'     estimates as a flat named vector.}
#'   \item{`settings`}{The arguments the fit was made with, including any
#'     user-supplied variances, aligned with `gaugings`. Enough to refit.}
#'   \item{`variance`}{For constructors with a `variance` argument: the
#'     variance scheme, with any estimated parameter filled in, such as the
#'     power estimated under [var_power()] (`fit$variance$exponent`).}
#' }
#'
#' Use [fitted()], [residuals()], [coef()] and [predict()] rather than
#' reaching into the fit for those quantities.
#'
#' Fits also carry elements that depend on how they are computed, and that
#' may change as the methods for modelling the scatter develop: `model`, the
#' underlying model object, such as from [stats::nls()]; `rse`, the residual
#' standard error (on the log scale for the log-scale fits); and, for
#' constructors with a `variance` argument, `weights_used`, the weights the
#' final model was fitted with, and `irls` (for iteratively reweighted least
#' squares): `NULL`, or under [var_prop()], which refits in rounds, a list
#' giving the number of rounds (`iterations`) and whether they `converged`.
#' Some constructors carry further elements, documented on their own help
#' pages.
#'
#' @name rating_curve
NULL


#' Print a rating curve fit
#'
#' Every fit carries `"rating_curve"` as its second class, so this one method
#' covers all of them.
#'
#' @param x A rating curve fit.
#' @param ... Ignored.
#' @return `x`, invisibly. Called for the printed output.
#' @export
print.rating_curve <- function(x, ...) {
  cl <- class(x)
  cl <- cl[cl != "rating_curve"]
  cl <- paste(cl, collapse = ", ")
  cat("Rating curve model.\n- Method:", cl, "\n")
  invisible(x)
}


#' Estimated parameters of a rating curve
#'
#' The contents of `object$curve_parameters` as a flat named numeric vector. For a
#' single-segment curve the names are those of `curve_parameters`; for a multi-segment
#' curve they carry the segment number, as in the model formula (`a1`, `b2`,
#' `k`, ...).
#'
#' @param object A rating curve fit.
#' @param ... Ignored.
#' @return A named numeric vector; empty for [rc_loess()].
#' @examples
#' coef(rc_powerlaw(discharge, stage, data = thompson))
#' @export
coef.rating_curve <- function(object, ...) {
  pars <- object$curve_parameters
  if (!length(pars)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  unlist(pars)
}


#' @rdname coef.rating_curve
#' @export
coef.rc_2seg_powerlaw <- function(object, ...) {
  stats::coef(object$model)
}


#' Fitted values and residuals of a rating curve
#'
#' The fitted discharge at each gauging, and the residuals: the gaugings'
#' departures from the curve.
#'
#' Residuals are calculated on the scale the model is fitted on: the log
#' scale for [rc_powerlaw_log()], and discharge, in cubic meters per second,
#' for every other fit. Nothing is transformed back. The `type` of residual
#' is one of:
#'
#' * `"difference"`: observed minus fitted, so `log(observed) - log(fitted)`
#'   for [rc_powerlaw_log()].
#' * `"scaled"`: the difference divided by the fit's estimated standard
#'   deviation of a gauging at that stage, on the same scale, as given by its
#'   variance scheme: for example `rse * fitted` under [var_prop()]. These
#'   are also known as Pearson residuals. If the variance scheme describes
#'   the scatter well, they have about the same spread at every stage, with a
#'   standard deviation near 1, so they are the ones to check for a trend in
#'   the scatter or for normality. The standard deviation is that of a
#'   gauging's scatter, not of the residual itself: it makes no adjustment
#'   for the uncertainty in the fitted curve.
#'
#' @param object A rating curve fit.
#' @param ... Must be empty.
#' @param type The type of residual: `"difference"` (the default) or
#'   `"scaled"`.
#' @return A numeric vector, one value per gauging in `object$gaugings`.
#' @examples
#' fit <- rc_powerlaw(discharge, stage, data = thompson, variance = "prop")
#' head(fitted(fit))
#' head(residuals(fit))
#'
#' # the residuals in standard deviations, to check the variance scheme
#' qqnorm(residuals(fit, type = "scaled"))
#' @export
fitted.rating_curve <- function(object, ...) {
  rlang::check_dots_empty()
  stats::predict(object, new_stage = object$gaugings$stage)$fit
}


#' @rdname fitted.rating_curve
#' @export
residuals.rating_curve <- function(
  object,
  ...,
  type = c("difference", "scaled")
) {
  rlang::check_dots_empty()
  type <- rlang::arg_match(type)
  observed <- object$gaugings$discharge
  fit <- stats::fitted(object)
  if (inherits(object, "rc_powerlaw_log")) {
    # fitted on the log scale, with the same standard deviation throughout
    residual <- log(observed) - log(fit)
    sd <- object$rse
  } else {
    # the weights are the reciprocals of the variances, up to the scale rse^2
    residual <- observed - fit
    sd <- object$rse / sqrt(object$weights_used)
  }
  if (type == "scaled") {
    residual <- residual / sd
  }
  residual
}
