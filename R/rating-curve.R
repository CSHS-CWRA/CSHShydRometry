# The structure shared by every fit, and the methods that rely on it.

#' Rating curve fits
#'
#' Every `rc_*()` constructor returns a list with class
#' `c("rc_<method>", "rating_curve")` and these elements:
#'
#' \describe{
#'   \item{`gaugings`}{A tibble of the gaugings used, with columns
#'     `discharge` and `stage`, after dropping any with a missing value.}
#'   \item{`pars`}{The estimated parameters of the curve, as a named list with
#'     one element per parameter type. For a multi-segment curve an element
#'     holds one value per segment (for `k`, per breakpoint); elements may
#'     differ in length where a parameter is fixed by how the segments join. Empty
#'     for [rc_loess()], which has no parameters. [coef()] returns the same
#'     estimates as a flat named vector.}
#'   \item{`settings`}{The arguments the fit was made with, including any
#'     user-supplied weights, aligned with `gaugings`. Enough to refit.}
#'   \item{`rse`}{The residual standard error (on the log scale for the
#'     log-scale fits).}
#'   \item{`model`}{The underlying model object, e.g. from [stats::nls()].}
#' }
#'
#' Constructors with a `wts` argument also return `weights`, the
#' weights the final model was fitted with, and `irls` (for iteratively
#' reweighted least squares): `NULL`, or under [wts_prop()], which refits in
#' rounds, a list giving the number of rounds (`iterations`) and whether they
#' `converged`; and `exponent`: under [wts_prop()], the power of the flow to
#' which the scatter is proportional, as given or as estimated, and `NULL`
#' otherwise. Some constructors carry
#' further elements, documented on their own help pages.
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
#' The contents of `object$pars` as a flat named numeric vector. For a
#' single-segment curve the names are those of `pars`; for a multi-segment
#' curve they carry the segment number, as in the model formula (`a1`, `b2`,
#' `k`, ...).
#'
#' @param object A rating curve fit.
#' @param ... Ignored.
#' @return A named numeric vector; empty for [rc_loess()].
#' @examples
#' coef(rc_power(discharge, stage, data = thompson))
#' @export
coef.rating_curve <- function(object, ...) {
  pars <- object$pars
  if (!length(pars)) {
    return(stats::setNames(numeric(0), character(0)))
  }
  unlist(pars)
}


#' @rdname coef.rating_curve
#' @export
coef.rc_2seg_power <- function(object, ...) {
  stats::coef(object$model)
}
