# Variance schemes: the models of scatter a fit can assume.

#' Variance schemes for rating-curve fits
#'
#' How the variance of the gaugings about the curve is modelled. Pass one of
#' these as the `variance` argument of [rc_powerlaw()], [rc_poly()],
#' [rc_loess()] or [rc_2seg_powerlaw()]. The fits weight each gauging by the
#' reciprocal of its variance.
#'
#' **Experimental:** the interface for modelling the scatter is likely to
#' grow, and may change.
#'
#' * `var_none()`: the variance is the same at every flow (ordinary least
#'   squares).
#' * `var_prop()`: the standard deviation is proportional to the flow, that
#'   is, a constant coefficient of variation. The variances, proportional to
#'   `fitted^2`, depend on the fit itself, so the fit is repeated in rounds:
#'   fit, recompute the weights from the fitted values, refit, each round
#'   starting from the previous round's estimates. The rounds stop once no
#'   fitted discharge changes by more than a fraction `tol` from one round to
#'   the next, or after `maxiter` rounds, with a warning.
#' * `var_power()`: the standard deviation is proportional to a power of the
#'   flow, with the power estimated along with the curve. This adds a
#'   parameter, so the fit is made by generalised least squares with
#'   [nlme::gnls()], starting from the fit under `var_prop()`. It is available
#'   in [rc_powerlaw()] only: that is a matter of implementation, as the other
#'   fitting functions are not made with [nlme::gnls()]. The estimate can be
#'   unstable with few gaugings; compare it with the fit under `var_prop()`.
#'   The estimated power is kept on the fit, in `fit$variance$exponent`. The
#'   limits from `predict()` treat it as known, so they leave out the
#'   uncertainty in the power.
#' * `var_spec()`: the relative variances of the gaugings are known,
#'   typically from their reported uncertainties, and given directly. Only
#'   their relative sizes matter: the fit estimates a common scale factor,
#'   the residual standard error, which is 1 when the gaugings scatter
#'   exactly as their variances say. A new gauging has no variance given, so
#'   prediction limits are returned as `NA`.
#'
#' The strings `"none"`, `"prop"` and `"power"` are shorthand for
#' `var_none()`, `var_prop()` and `var_power()` with their defaults.
#'
#' @param values The variances, one per gauging, as a numeric vector: for
#'   example `sd^2`, where `sd` is the standard uncertainty of each discharge.
#'   They are evaluated where `var_spec()` is called, as an ordinary argument,
#'   not looked up in the fit's `data`; write `sauze$uncertainty_sd^2`, not
#'   `uncertainty_sd^2`. The variances of gaugings that the fit drops, for a
#'   missing stage or discharge, are dropped with them.
#' @param tol Convergence tolerance: the reweighting stops once no fitted
#'   discharge changes by more than this fraction from one round to the next.
#' @param maxiter Maximum number of reweighting rounds.
#' @return An object of class `"rc_variance"`.
#' @examples
#' rc_powerlaw(discharge, stage, data = thompson, variance = var_prop())
#'
#' # the same, with the defaults
#' rc_powerlaw(discharge, stage, data = thompson, variance = "prop")
#'
#' # the power of the flow estimated too
#' fit <- rc_powerlaw(discharge, stage, data = thompson, variance = var_power())
#' fit$variance$exponent
#'
#' # variances from each gauging's reported uncertainty, leaving out the one
#' # reported as 0.0226%, probably a fraction entered as a percentage (see
#' # ?thompson): it would carry over 99.9% of the total weight
#' d <- thompson[!is.na(thompson$uncertainty_pct) & thompson$uncertainty_pct > 1, ]
#' d$uncertainty_sd <- d$uncertainty_pct / 100 * d$discharge / 2
#' rc_powerlaw(discharge, stage, data = d, variance = var_spec(d$uncertainty_sd^2))
#' @name variance
NULL


#' @rdname variance
#' @export
var_none <- function() {
  new_variance("none")
}


#' @rdname variance
#' @export
var_prop <- function(tol = 1e-6, maxiter = 100) {
  checkmate::assert_number(tol, lower = 0)
  checkmate::assert_count(maxiter, positive = TRUE)
  # a known variance function: no parameter to estimate, so the fitting
  # functions handle it themselves, by reweighting in rounds
  var_nlme(nlme::varPower(fixed = 1), "prop", tol = tol, maxiter = maxiter)
}


#' @rdname variance
#' @export
var_power <- function() {
  # a variance function with a parameter to estimate, by nlme::gnls()
  var_nlme(nlme::varPower(), "power")
}


#' @rdname variance
#' @export
var_spec <- function(values) {
  checkmate::assert_numeric(values, min.len = 1L, lower = 0)
  if (any(values == 0, na.rm = TRUE)) {
    stop("The values of `var_spec()` must be positive.", call. = FALSE)
  }
  new_variance("spec", values = values)
}


#' @param x An `"rc_variance"` object.
#' @param ... Ignored.
#' @rdname variance
#' @export
print.rc_variance <- function(x, ...) {
  desc <- switch(
    x$type,
    none = "none (the same at every flow)",
    prop = sprintf(
      "standard deviation proportional to the flow (tol = %g, maxiter = %d)",
      x$tol,
      as.integer(x$maxiter)
    ),
    power = if (is.null(x$exponent)) {
      "standard deviation proportional to a power of the flow, the power to be estimated"
    } else {
      sprintf(
        "standard deviation proportional to the flow to the power %.3g (estimated)",
        x$exponent
      )
    },
    spec = sprintf("specified, for %d gaugings", length(x$values))
  )
  cat("Variance:", desc, "\n")
  invisible(x)
}


#' Construct a variance scheme
#'
#' @param type `"none"`, `"prop"`, `"power"` or `"spec"`.
#' @param ... Its settings.
#' @return An object of class `c("rc_var_<type>", "rc_variance")`.
#' @noRd
new_variance <- function(type, ...) {
  structure(
    list(type = type, ...),
    class = c(paste0("rc_var_", type), "rc_variance")
  )
}


#' A variance scheme described by an nlme variance function
#'
#' The variance as a function of the fitted values, in the terms of
#' [nlme::varFunc()], for `var_power()`, whose power [nlme::gnls()]
#' estimates; `var_prop()` is built on it too, for consistency. This is
#' internal, and not a commitment to nlme: if the estimation moves in-house,
#' it can go. A variance function with no free parameters,
#' like `var_prop()`'s, need not be fitted by [nlme::gnls()]: the fitting
#' functions reweight in rounds instead.
#'
#' @param varfunc An nlme variance function, such as [nlme::varPower()].
#' @param type The scheme's name, for printing and dispatch.
#' @param ... Further settings, such as those of the reweighting.
#' @return An `"rc_variance"` object.
#' @noRd
var_nlme <- function(varfunc, type, ...) {
  checkmate::assert_class(varfunc, "varFunc")
  new_variance(type, varfunc = varfunc, ...)
}


#' Resolve the `variance` argument of a fitting function
#'
#' Accepts an `"rc_variance"` object or the shorthand `"none"`, `"prop"` or
#' `"power"`, and, for specified variances, drops those of gaugings with a
#' missing stage or discharge, so that the variances stay aligned with the
#' gaugings the fit keeps.
#'
#' @param variance The argument as given.
#' @param keep Logical vector: which gaugings are kept.
#' @param fitter The fitting function, for error messages.
#' @param power_ok Whether the fitting function supports [var_power()].
#' @return An `"rc_variance"` object.
#' @noRd
resolve_variance <- function(variance, keep, fitter, power_ok = FALSE) {
  if (is.character(variance)) {
    variance <- switch(
      rlang::arg_match0(
        variance,
        c("none", "prop", "power", "spec"),
        arg_nm = "variance"
      ),
      none = var_none(),
      prop = var_prop(),
      power = var_power(),
      spec = stop(
        "`variance = \"spec\"` needs the variances: use `var_spec(values)`.",
        call. = FALSE
      )
    )
  }
  if (!inherits(variance, "rc_variance")) {
    stop(
      "`variance` must be \"none\", \"prop\", \"power\", or made by ",
      "`var_none()`, `var_prop()`, `var_power()` or `var_spec()`.",
      call. = FALSE
    )
  }
  if (variance$type == "power" && !power_ok) {
    stop(
      "`var_power()` is available in `rc_powerlaw()` only, not in `", fitter,
      "()`: its power is estimated with `nlme::gnls()`, which `", fitter,
      "()` is not fitted with. Use `var_prop()` for a standard deviation ",
      "proportional to the flow.",
      call. = FALSE
    )
  }
  if (variance$type == "spec") {
    checkmate::assert_numeric(
      variance$values,
      len = length(keep),
      .var.name = "the values of var_spec()"
    )
    variance$values <- variance$values[keep]
  }
  variance
}
