# Weighting schemes: the error models a fit can assume.

#' Weighting schemes for rating-curve fits
#'
#' How the scatter of the gaugings about the curve is modelled. Pass one of
#' these as the `wts` argument of [rc_power()], [rc_poly()], [rc_loess()] or
#' [rc_2seg_power()].
#'
#' * `wts_none()`: the scatter is the same at every flow (ordinary least
#'   squares).
#' * `wts_prop()`: the scatter is proportional to the flow, that is, a
#'   constant coefficient of variation. The weights, `1 / fitted^2`, depend on
#'   the fit itself, so the fit is repeated in rounds: fit, recompute the
#'   weights from the fitted values, refit, each round starting from the
#'   previous round's estimates. The rounds stop once no fitted discharge
#'   changes by more than a fraction `tol` from one round to the next, or
#'   after `maxiter` rounds, with a warning.
#' * `wts_power()`: the scatter is proportional to a power of the flow, with
#'   the power estimated along with the curve. This adds a parameter, so the
#'   fit is made by generalised least squares with [nlme::gnls()], starting
#'   from the fit under `wts_prop()`. It is available in [rc_power()] only:
#'   that is a matter of implementation, as the other fitting functions are not
#'   made with [nlme::gnls()]. The estimate can be unstable with few gaugings;
#'   compare it with the fit under `wts_prop()`. The estimated power is kept on
#'   the fit, in `fit$wts$exponent`.
#' * `wts_spec()`: the scatter of each gauging is known, typically from its
#'   reported uncertainty, and given as weights: the reciprocal of each
#'   gauging's variance. A new gauging's scatter is then not estimated, so
#'   prediction limits are returned as `NA`.
#'
#' The strings `"none"`, `"prop"` and `"power"` are shorthand for
#' `wts_none()`, `wts_prop()` and `wts_power()` with their defaults.
#'
#' @param values <[`data-masking`][rlang::args_data_masking]> The weights, one
#'   per gauging: a vector, or an expression evaluated in the fit's `data`,
#'   such as `1 / uncertainty_sd^2`.
#' @param tol Convergence tolerance: the reweighting stops once no fitted
#'   discharge changes by more than this fraction from one round to the next.
#' @param maxiter Maximum number of reweighting rounds.
#' @return An object of class `"rc_wts"`.
#' @examples
#' rc_power(discharge, stage, data = thompson, wts = wts_prop())
#'
#' # the same, with the defaults
#' rc_power(discharge, stage, data = thompson, wts = "prop")
#'
#' # the power of the flow estimated too
#' fit <- rc_power(discharge, stage, data = thompson, wts = wts_power())
#' fit$wts$exponent
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
wts_prop <- function(tol = 1e-6, maxiter = 100) {
  checkmate::assert_number(tol, lower = 0)
  checkmate::assert_count(maxiter, positive = TRUE)
  # a known variance function: no parameter to estimate, so the fitting
  # functions handle it themselves, by reweighting in rounds
  wts_nlme(nlme::varPower(fixed = 1), "prop", tol = tol, maxiter = maxiter)
}


#' @rdname wts
#' @export
wts_power <- function() {
  # a variance function with a parameter to estimate, by nlme::gnls()
  wts_nlme(nlme::varPower(), "power")
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
    prop = sprintf(
      "proportional to the flow (tol = %g, maxiter = %d)",
      x$tol,
      as.integer(x$maxiter)
    ),
    power = if (is.null(x$exponent)) {
      "proportional to a power of the flow, the power to be estimated"
    } else {
      sprintf(
        "proportional to the flow to the power %.3g (estimated)",
        x$exponent
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
#' @param type `"none"`, `"prop"`, `"power"` or `"spec"`.
#' @param ... Its settings.
#' @return An object of class `c("rc_wts_<type>", "rc_wts")`.
#' @noRd
new_wts <- function(type, ...) {
  structure(
    list(type = type, ...),
    class = c(paste0("rc_wts_", type), "rc_wts")
  )
}


#' A weighting scheme described by an nlme variance function
#'
#' The scatter as a function of the fitted values, in the terms of
#' [nlme::varFunc()]. `wts_prop()` and `wts_power()` are built on this, so
#' that a public version could later accept any of nlme's variance functions,
#' fitted by [nlme::gnls()]. A variance function with no free parameters,
#' like `wts_prop()`'s, need not be fitted by [nlme::gnls()]: the fitting
#' functions reweight in rounds instead.
#'
#' @param variance An nlme variance function, such as [nlme::varPower()].
#' @param type The scheme's name, for printing and dispatch.
#' @param ... Further settings, such as those of the reweighting.
#' @return An `"rc_wts"` object.
#' @noRd
wts_nlme <- function(variance, type, ...) {
  checkmate::assert_class(variance, "varFunc")
  new_wts(type, variance = variance, ...)
}


#' Resolve the `wts` argument of a fitting function
#'
#' Accepts an `"rc_wts"` object or the shorthand `"none"`, `"prop"` or
#' `"power"`, and, for specified weights, evaluates the values in `data` and
#' drops those of gaugings with a missing stage or discharge. The result
#' holds the values themselves, so it can be stored on the fit and reused to
#' refit.
#'
#' @param wts The argument as given.
#' @param data The fit's `data`, or `NULL`.
#' @param keep Logical vector: which gaugings are kept.
#' @param fitter The fitting function, for error messages.
#' @param power_ok Whether the fitting function supports [wts_power()].
#' @return An `"rc_wts"` object.
#' @noRd
resolve_wts <- function(wts, data, keep, fitter, power_ok = FALSE) {
  if (is.character(wts)) {
    wts <- switch(
      rlang::arg_match0(
        wts,
        c("none", "prop", "power", "spec"),
        arg_nm = "wts"
      ),
      none = wts_none(),
      prop = wts_prop(),
      power = wts_power(),
      spec = stop(
        "`wts = \"spec\"` needs the weights: use `wts_spec(values)`.",
        call. = FALSE
      )
    )
  }
  if (!inherits(wts, "rc_wts")) {
    stop(
      "`wts` must be \"none\", \"prop\", \"power\", or made by `wts_none()`, ",
      "`wts_prop()`, `wts_power()` or `wts_spec()`.",
      call. = FALSE
    )
  }
  if (wts$type == "power" && !power_ok) {
    stop(
      "`wts_power()` is available in `rc_power()` only, not in `", fitter,
      "()`: its power is estimated with `nlme::gnls()`, which `", fitter,
      "()` is not fitted with. Use `wts_prop()` for scatter proportional ",
      "to the flow.",
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
