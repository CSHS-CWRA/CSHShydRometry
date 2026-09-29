# Curve evaluation and the parameter scale used to
# propagate uncertainty.

#' Evaluate the two-segment curve
#'
#' The mean function itself, in ordinary (unlogged) parameters. Every caller
#' arrives here eventually, whichever parameter scale it works on: the scales
#' differ only in how the parameters are recovered, not in what the curve is,
#' so this is the one place the model is written down for evaluation.
#'
#' Vectorised over parameter vectors rather than over stages: pass a scalar
#' `hh` and vectors of parameters (one entry per draw) to get one discharge
#' per draw. `ifelse()` recycles the scalar `hh` against them.
#'
#' @param hh Stage at which to evaluate, normally a single value.
#' @param a1,b1,c1 Lower-segment coefficient, exponent and stage offset.
#' @param k Breakpoint stage.
#' @param b2 Upper-segment exponent.
#' @param c2 Upper-segment offset. Piecewise only; the upper coefficient is
#'   not free there, being fixed by continuity at `h = k`.
#' @param a2 Upper-segment coefficient. Compound only.
#' @param cfg `"piecewise"` or `"compound"`.
#' @return Numeric vector, as long as the parameter vectors.
#' @keywords internal
rc_qeval <- function(hh, a1, b1, c1, k, b2, c2 = NULL, a2 = NULL, cfg) {
  if (cfg == "piecewise") {
    ifelse(
      hh < k,
      a1 * (hh - c1)^b1,
      (a1 * (k - c1)^b1 / (k - c2)^b2) * (hh - c2)^b2
    )
  } else {
    ifelse(hh < k, a1 * (hh - c1)^b1, a1 * (k - c1)^b1 + a2 * (hh - k)^b2)
  }
}


#' Parameter vector, covariance and curve evaluator on a chosen scale
#'
#' Rating-curve parameters are constrained: `a` and `b` are positive, the offsets
#' must lie below the stages at which the curve is evaluated, and `c2 < k`. A
#' normal distribution on the fitted parameters respects none of this, so drawing
#' from it can land where the curve is undefined. Re-expressing the model in
#' terms of the gaps between parameters, on a log scale, removes the problem:
#' every point of R^p then maps to a valid curve.
#'
#' `space = "original"` returns the fit as-is. `space = "log"` refits the same
#' model in the coordinates
#' \code{log a1, log b1, log(hmin - c1), log(k - c2) [or log a2], log b2,
#' log(k - hmin)}, where `hmin` is the smallest observed stage. The fitted curve
#' is unchanged -- it is the same model in different coordinates -- but the
#' normal approximation is then formed on a scale where it cannot escape the
#' parameter space.
#'
#' @param object An `rc_nls_2seg` fit.
#' @param space `"original"` or `"log"`.
#' @return A list with `th` (parameters), `Sg` (covariance), `ik` (name of the
#'   breakpoint parameter on this scale) and `qmat(hh, P)`, which evaluates the
#'   curve at stage `hh` for a matrix `P` of parameter vectors (one per row).
#' @keywords internal
rc_param_space <- function(object, space = c("original", "log")) {
  space <- rlang::arg_match(space)
  cfg <- object$configuration
  hmin <- min(object$qh_obs$h)
  if (space == "original") {
    th <- stats::coef(object$model)
    Sg <- stats::vcov(object$model)
    # On this scale the columns of P are already the model parameters. The
    # config is fixed once here rather than tested on every
    # call, since sim_limits_2seg() evaluates it once per stage
    # over many draws.
    qmat <- if (cfg == "piecewise") {
      function(hh, P) {
        rc_qeval(
          hh, P[, "a1"], P[, "b1"], P[, "c1"], P[, "k"], P[, "b2"],
          c2 = P[, "c2"], cfg = "piecewise"
        )
      }
    } else {
      function(hh, P) {
        rc_qeval(
          hh, P[, "a1"], P[, "b1"], P[, "c1"], P[, "k"], P[, "b2"],
          a2 = P[, "a2"], cfg = "compound"
        )
      }
    }
    return(list(th = th, Sg = Sg, ik = "k", qmat = qmat))
  }
  # --- log scale: refit the same model in gap coordinates ---
  qv <- object$qh_obs$q
  hv <- object$qh_obs$h
  th0 <- stats::coef(object$model)
  w <- object$weights
  if (cfg == "piecewise") {
    st <- list(
      la1 = log(th0[["a1"]]),
      lb1 = log(th0[["b1"]]),
      lgam = log(hmin - th0[["c1"]]),
      lb2 = log(th0[["b2"]]),
      ldel = log(th0[["k"]] - th0[["c2"]]),
      lkap = log(th0[["k"]] - hmin)
    )
    form <- q ~ ifelse(
      h < hmin + exp(lkap),
      exp(la1) * (h - (hmin - exp(lgam)))^exp(lb1),
      (exp(la1) *
        (hmin + exp(lkap) - (hmin - exp(lgam)))^exp(lb1) /
        exp(ldel)^exp(lb2)) *
        (h - (hmin + exp(lkap) - exp(ldel)))^exp(lb2)
    )
  } else {
    st <- list(
      la1 = log(th0[["a1"]]),
      lb1 = log(th0[["b1"]]),
      lgam = log(hmin - th0[["c1"]]),
      la2 = log(th0[["a2"]]),
      lb2 = log(th0[["b2"]]),
      lkap = log(th0[["k"]] - hmin)
    )
    form <- q ~ ifelse(
      h < hmin + exp(lkap),
      exp(la1) * (h - (hmin - exp(lgam)))^exp(lb1),
      exp(la1) *
        (hmin + exp(lkap) - (hmin - exp(lgam)))^exp(lb1) +
        exp(la2) * (h - (hmin + exp(lkap)))^exp(lb2)
    )
  }
  m <- nls2::nls2(
    form,
    data = data.frame(q = qv, h = hv),
    weights = w,
    start = st,
    algorithm = "port",
    control = list(tol = 1e-6, maxiter = 500)
  )
  # Here the columns of P are gap coordinates, so undo the reparameterisation
  # first; the curve itself is the same one rc_qeval() always evaluates.
  qmat <- if (cfg == "piecewise") {
    function(hh, P) {
      k <- hmin + exp(P[, "lkap"])
      rc_qeval(
        hh, exp(P[, "la1"]), exp(P[, "lb1"]), hmin - exp(P[, "lgam"]), k,
        exp(P[, "lb2"]), c2 = k - exp(P[, "ldel"]), cfg = "piecewise"
      )
    }
  } else {
    function(hh, P) {
      k <- hmin + exp(P[, "lkap"])
      rc_qeval(
        hh, exp(P[, "la1"]), exp(P[, "lb1"]), hmin - exp(P[, "lgam"]), k,
        exp(P[, "lb2"]), a2 = exp(P[, "la2"]), cfg = "compound"
      )
    }
  }
  list(th = stats::coef(m), Sg = stats::vcov(m), ik = "lkap", qmat = qmat)
}


#' Evaluate the curve for a single parameter vector
#'
#' Convenience wrapper on `sp$qmat()`, which expects a matrix of parameter
#' vectors: this wraps one named vector as a single row. Handy when writing
#' the calculation out by hand, as `derivation/reproduce_results.R` does.
#' Prefer passing a whole matrix where you can -- the interval methods
#' evaluate thousands of parameter vectors, and one call for all of them is
#' far cheaper than one call each.
#'
#' @param sp A list from [rc_param_space()].
#' @param hh Stage at which to evaluate.
#' @param p A named parameter vector, on the scale `sp` was built for.
#' @return A single discharge value.
#' @keywords internal
rc_qpoint <- function(sp, hh, p) {
  as.numeric(sp$qmat(hh, matrix(p, nrow = 1, dimnames = list(NULL, names(p)))))
}
