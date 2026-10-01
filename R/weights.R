# Weighting schemes: the error models a fit can assume.

#' Weighting schemes for rating-curve fits
#'
#' How the scatter of the gaugings about the curve is modelled. Pass one of
#' these as the `wts` argument of [rc_power()], [rc_poly()], [rc_loess()] or
#' [rc_2seg_power()].
#'
#' * `wts_none()`: the scatter is the same at every flow (ordinary least
#'   squares).
#' * `wts_prop()`: the scatter is proportional to the flow raised to
#'   `exponent`; with the default `exponent = 1`, a constant coefficient of
#'   variation. The weights, `1 / fitted^(2 * exponent)`, depend on the fit
#'   itself, so the fit is repeated in rounds: fit, recompute the weights from
#'   the fitted values, refit, each round starting from the previous round's
#'   estimates. The rounds stop once no fitted discharge changes by more than a
#'   fraction `tol` from one round to the next, or after `maxiter` rounds, with
#'   a warning.
#'
#'   With `exponent = NULL`, the exponent is estimated along with the curve,
#'   by generalised least squares with [nlme::gnls()], starting from the fit
#'   with `exponent = 1`. This is available in [rc_power()] only. The estimate
#'   can be unstable with few gaugings: check it against fits with the
#'   exponent fixed.
#' * `wts_spec()`: the scatter of each gauging is known, typically from its
#'   reported uncertainty, and given as weights: the reciprocal of each
#'   gauging's variance. A new gauging's scatter is then not estimated, so
#'   prediction limits are returned as `NA`.
#'
#' The strings `"none"` and `"prop"` are shorthand for `wts_none()` and
#' `wts_prop()` with its defaults.
#'
#' @param values <[`data-masking`][rlang::args_data_masking]> The weights, one
#'   per gauging: a vector, or an expression evaluated in the fit's `data`,
#'   such as `1 / uncertainty_sd^2`.
#' @param tol Convergence tolerance: the reweighting stops once no fitted
#'   discharge changes by more than this fraction from one round to the next.
#' @param maxiter Maximum number of reweighting rounds.
#' @param exponent The power of the flow to which the scatter is
#'   proportional: a number (1 by default), or `NULL` to estimate it.
#' @return An object of class `"rc_wts"`.
#' @examples
#' rc_power(discharge, stage, data = thompson, wts = wts_prop())
#'
#' # the same, with the defaults
#' rc_power(discharge, stage, data = thompson, wts = "prop")
#'
#' # scatter growing faster than the flow, or with the exponent estimated
#' rc_power(discharge, stage, data = thompson, wts = wts_prop(exponent = 1.5))
#' rc_power(discharge, stage, data = thompson, wts = wts_prop(exponent = NULL))
#'
#' # weights from each gauging's reported uncertainty
#' d <- thompson[!is.na(thompson$uncertainty_pct), ]
#' d$uncertainty_sd <- d$uncertainty_pct / 100 * d$discharge / 2
#' rc_power(discharge, stage, data = d, wts = wts_spec(1 / uncertainty_sd^2))
#' @name wts
NULL


#' @rdname wts
#' @export
wts_none <- function() {
  new_wts("none")
}


#' @rdname wts
#' @export
wts_prop <- function(tol = 1e-6, maxiter = 100, exponent = 1) {
  checkmate::assert_number(tol, lower = 0)
  checkmate::assert_count(maxiter, positive = TRUE)
  checkmate::assert_number(exponent, null.ok = TRUE, finite = TRUE)
  new_wts("prop", tol = tol, maxiter = maxiter, exponent = exponent)
}


#' @rdname wts
#' @export
wts_spec <- function(values) {
  new_wts("spec", values = rlang::enquo(values))
}


#' @param x An `"rc_wts"` object.
#' @param ... Ignored.
#' @rdname wts
#' @export
print.rc_wts <- function(x, ...) {
  desc <- switch(
    x$type,
    none = "none (the same scatter at every flow)",
    prop = if (is.null(x$exponent)) {
      "proportional to a power of the flow, the exponent estimated"
    } else if (x$exponent == 1) {
      sprintf(
        "proportional to the flow (tol = %g, maxiter = %d)",
        x$tol,
        as.integer(x$maxiter)
      )
    } else {
      sprintf(
        "proportional to the flow to the power %g (tol = %g, maxiter = %d)",
        x$exponent,
        x$tol,
        as.integer(x$maxiter)
      )
    },
    spec = if (rlang::is_quosure(x$values)) {
      paste("specified:", rlang::as_label(x$values))
    } else {
      sprintf("specified, for %d gaugings", length(x$values))
    }
  )
  cat("Weighting:", desc, "\n")
  invisible(x)
}


#' Construct a weighting scheme
#'
#' @param type `"none"`, `"prop"` or `"spec"`.
#' @param ... Its settings.
#' @return An object of class `c("rc_wts_<type>", "rc_wts")`.
#' @noRd
new_wts <- function(type, ...) {
  structure(
    list(type = type, ...),
    class = c(paste0("rc_wts_", type), "rc_wts")
  )
}


#' Resolve the `wts` argument of a fitting function
#'
#' Accepts an `"rc_wts"` object or the shorthand `"none"` or `"prop"`, and,
#' for specified weights, evaluates the values in `data` and drops those of
#' gaugings with a missing stage or discharge. The result holds the values
#' themselves, so it can be stored on the fit and reused to refit.
#'
#' @param wts The argument as given.
#' @param data The fit's `data`, or `NULL`.
#' @param keep Logical vector: which gaugings are kept.
#' @param estimate_exponent Whether the fit can estimate the exponent of
#'   [wts_prop()].
#' @return An `"rc_wts"` object.
#' @noRd
resolve_wts <- function(wts, data, keep, estimate_exponent = FALSE) {
  if (is.character(wts)) {
    wts <- switch(
      rlang::arg_match0(wts, c("none", "prop", "spec"), arg_nm = "wts"),
      none = wts_none(),
      prop = wts_prop(),
      spec = stop(
        "`wts = \"spec\"` needs the weights: use `wts_spec(values)`.",
        call. = FALSE
      )
    )
  }
  if (!inherits(wts, "rc_wts")) {
    stop(
      "`wts` must be \"none\", \"prop\", or made by `wts_none()`, ",
      "`wts_prop()` or `wts_spec()`.",
      call. = FALSE
    )
  }
  if (wts$type == "prop" && is.null(wts$exponent) && !estimate_exponent) {
    stop(
      "Estimating the exponent of `wts_prop()` is available in `rc_power()` ",
      "only; give a value for `exponent` here.",
      call. = FALSE
    )
  }
  if (wts$type == "spec") {
    values <- rlang::eval_tidy(wts$values, data)
    checkmate::assert_numeric(values, len = length(keep), .var.name = "wts")
    wts$values <- values[keep]
  }
  wts
}
