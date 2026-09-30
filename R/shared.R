# Helpers shared by the one- and two-segment models.

#' Confidence and prediction limits for an NLS fit with proportional weights
#'
#' `nlspw` is short for "NLS, proportional weights": these are the limits
#' used when `wts_code = "prop"`.
#'
#' Builds pointwise interval limits for the `wts_code = "prop"` case, where the
#' error standard deviation is taken to be proportional to the mean discharge
#' (constant coefficient of variation). Standard errors of the fitted mean come
#' from the delta method via [investr::predFit()]; the interval half-width uses
#' a Student-t multiplier.
#'
#' @param mod A fitted `nls` model object.
#' @param type `"confidence"` (uncertainty in the mean curve only) or
#'   `"prediction"` (adds the scatter of an individual observation).
#' @param level Coverage probability, e.g. `0.95`.
#' @param hpred Numeric vector of stage values at which to evaluate the limits.
#' @return A data frame with columns `fit`, `lwr`, `upr`.
#' @details For a confidence interval the half-width is `t * se(fit)`. For a
#'   prediction interval the observation variance `residual.scale^2 / w` is
#'   added, where the prediction-point weight `w = 1/fit^2` reintroduces the
#'   proportional-error assumption (so the added variance grows as `fit^2`).
#' @keywords internal
nlspw_limits <- function(
  mod,
  hpred,
  type = c("confidence", "prediction"),
  level = 0.95
) {
  type <- rlang::arg_match(type)

  ## compute fitted values, se and residual scale for hpred values
  nlsw_predfit <- investr::predFit(
    mod,
    newdata = data.frame(h = hpred),
    se.fit = TRUE
  )
  nlsw_se <- nlsw_predfit$se.fit
  nlsw_resscale <- nlsw_predfit$residual.scale
  yp_nlsw <- nlsw_predfit$fit
  nlsw_df <- nlsw_predfit$df

  ## compute weights for predicted values
  wtsp <- 1 / yp_nlsw^2

  ## compute limits for confidence and prediction intervals
  if (type == "confidence") {
    sp <- abs(nlsw_se)
  } else {
    ## prediction limits
    sp <- sqrt(nlsw_se^2 + (nlsw_resscale^2 / wtsp))
  }
  tc <- stats::qt(p = 0.5 + 0.5 * level, df = nlsw_df)
  nlsw_lims <- data.frame(
    fit = yp_nlsw,
    lwr = yp_nlsw - tc * sp,
    upr = yp_nlsw + tc * sp
  )
  nlsw_lims
}


#' Default stage grid for prediction
#'
#' The grid used when `hpred` is not supplied: 1000 points spanning the
#' observed stage range. Factored out so that every function taking `hpred`
#' treats a missing one the same way -- extrapolating beyond the gaugings is a
#' decision for the caller to make deliberately, not a default.
#'
#' @param object A fitted rating curve carrying `qh_obs`.
#' @param n Number of grid points.
#' @return Numeric vector of stage values.
#' @keywords internal
rc_hpred_grid <- function(object, n = 1000) {
  h <- object[["qh_obs"]][["h"]]
  checkmate::assert_numeric(h, min.len = 1L, any.missing = FALSE)
  seq(min(h), max(h), length.out = n)
}


#' Drop gaugings with a missing stage or discharge
#'
#' @param q,h Discharge and stage vectors of equal length.
#' @return A data frame with columns `q` and `h`, holding the complete cases.
#' @keywords internal
rc_complete <- function(q, h) {
  keep <- stats::complete.cases(q, h)
  data.frame(q = q[keep], h = h[keep])
}



#' Fit with proportional weights by iterative reweighting
#'
#' Under `wts_code = "prop"` the error standard deviation is proportional to
#' the mean discharge, so the weights `1 / fitted^2` depend on the fit itself.
#' The fit is therefore repeated in rounds: fit with the current weights,
#' recompute the weights from the new fitted values, refit. Each round starts
#' from the previous round's estimates. The rounds stop once no fitted
#' discharge changes by more than a fraction `wts_tol` between successive
#' rounds, or after `wts_maxiter` rounds, with a warning.
#'
#' @param fit_fun Function of `(wts, start)` returning a fitted model with
#'   `predict()` and `coef()` methods.
#' @param yp Initial fitted discharges, from which the first weights are
#'   computed. Must be positive.
#' @param start Named list of starting values for the first round.
#' @param wts_tol,wts_maxiter Convergence tolerance and maximum number of
#'   rounds.
#' @return A list with the final `model`, the `weights` it was fitted with,
#'   and `irls`: a list of the number of `iterations` (rounds) and whether
#'   the rounds `converged`.
#' @keywords internal
rc_irls <- function(fit_fun, yp, start, wts_tol, wts_maxiter) {
  converged <- FALSE
  for (i in seq_len(wts_maxiter)) {
    wts <- 1 / yp^2
    mod <- fit_fun(wts, start)
    yp_new <- as.numeric(stats::predict(mod))
    change <- max(abs(yp_new - yp) / abs(yp))
    start <- as.list(stats::coef(mod))
    yp <- yp_new
    if (change < wts_tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    warning(
      sprintf(
        paste(
          "Proportional weights did not converge in %d rounds (the fitted",
          "discharges still changed by up to %.2g%% in the last); the fit",
          "may not be reliable. Consider increasing `wts_maxiter`."
        ),
        wts_maxiter,
        100 * change
      ),
      call. = FALSE
    )
  }
  list(
    model = mod,
    weights = wts,
    irls = list(iterations = i, converged = converged)
  )
}
