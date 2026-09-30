# Helpers shared by the one- and two-segment models.

#' Confidence and prediction limits for an NLS fit with proportional weights
#'
#' `nlspw` is short for "NLS, proportional weights": these are the limits
#' used under [wts_prop()].
#'
#' Builds pointwise interval limits for the [wts_prop()] case, where the
#' error standard deviation is taken to be proportional to the mean discharge
#' (constant coefficient of variation). Standard errors of the fitted mean come
#' from the delta method via [investr::predFit()]; the interval half-width uses
#' a Student-t multiplier.
#'
#' @param mod A fitted `nls` model object.
#' @param type `"confidence"` (uncertainty in the mean curve only) or
#'   `"prediction"` (adds the scatter of an individual observation).
#' @param level Coverage probability, e.g. `0.95`.
#' @param stage Numeric vector of stages at which to evaluate the limits.
#' @return A data frame with columns `fit`, `lwr`, `upr`.
#' @details For a confidence interval the half-width is `t * se(fit)`. For a
#'   prediction interval the observation variance `residual.scale^2 / w` is
#'   added, where the prediction-point weight `w = 1/fit^2` reintroduces the
#'   proportional-error assumption (so the added variance grows as `fit^2`).
#' @keywords internal
nlspw_limits <- function(
  mod,
  stage,
  type = c("confidence", "prediction"),
  level = 0.95
) {
  type <- rlang::arg_match(type)

  ## compute fitted values, se and residual scale for stage values
  nlsw_predfit <- investr::predFit(
    mod,
    newdata = data.frame(stage = stage),
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
#' The grid used when `stage` is not supplied: 1000 points spanning the
#' observed stage range. Factored out so that every function taking `stage`
#' treats a missing one the same way -- extrapolating beyond the gaugings is a
#' decision for the caller to make deliberately, not a default.
#'
#' @param object A fitted rating curve carrying `gaugings`.
#' @param n Number of grid points.
#' @return Numeric vector of stage values.
#' @keywords internal
rc_stage_grid <- function(object, n = 1000) {
  stage <- object[["gaugings"]][["stage"]]
  checkmate::assert_numeric(stage, min.len = 1L, any.missing = FALSE)
  seq(min(stage), max(stage), length.out = n)
}


#' Drop gaugings with a missing stage or discharge
#'
#' @param discharge,stage Vectors of equal length.
#' @return A data frame with columns `discharge` and `stage`, holding the
#'   complete cases.
#' @keywords internal
rc_complete <- function(discharge, stage) {
  keep <- stats::complete.cases(discharge, stage)
  data.frame(discharge = discharge[keep], stage = stage[keep])
}



#' Fit with proportional weights by iterative reweighting
#'
#' Under [wts_prop()] the error standard deviation is proportional to
#' the mean discharge, so the weights `1 / fitted^2` depend on the fit itself.
#' The fit is therefore repeated in rounds: fit with the current weights,
#' recompute the weights from the new fitted values, refit. Each round starts
#' from the previous round's estimates. The rounds stop once no fitted
#' discharge changes by more than a fraction `tol` between successive
#' rounds, or after `maxiter` rounds, with a warning.
#'
#' @param fit_fun Function of `(wts, start)` returning a fitted model with a
#'   `predict()` method, and a `coef()` method unless `start` is `NULL`.
#' @param yp Initial fitted discharges, from which the first weights are
#'   computed. Must be positive.
#' @param start Named list of starting values for the first round, or
#'   `NULL` for a model, such as loess, that needs none.
#' @param tol,maxiter Convergence tolerance and maximum number of rounds.
#' @return A list with the final `model`, the `weights` it was fitted with,
#'   and `irls`: a list of the number of `iterations` (rounds) and whether
#'   the rounds `converged`.
#' @keywords internal
rc_irls <- function(fit_fun, yp, start, tol, maxiter) {
  converged <- FALSE
  for (i in seq_len(maxiter)) {
    wts <- 1 / yp^2
    mod <- fit_fun(wts, start)
    yp_new <- as.numeric(stats::predict(mod))
    change <- max(abs(yp_new - yp) / abs(yp))
    if (!is.null(start)) {
      start <- as.list(stats::coef(mod))
    }
    yp <- yp_new
    if (change < tol) {
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
          "may not be reliable. Consider increasing `maxiter` in `wts_prop()`."
        ),
        maxiter,
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


#' Profile log-likelihood of a fit under its error model
#'
#' The normal log-likelihood, up to a constant, with the error variance of
#' each observation proportional to `1 / w` and the common scale profiled
#' out. Used to choose between fits of the same model from different starting
#' values: with fixed weights it orders fits exactly as the weighted residual
#' sum of squares does, and it stays comparable when, as under proportional
#' weights, the weights themselves depend on the fit.
#'
#' @param discharge Observed discharges.
#' @param mu Fitted discharges.
#' @param w Weights, the reciprocal of each observation's relative variance.
#' @return A single number.
#' @keywords internal
rc_loglik <- function(discharge, mu, w) {
  n <- length(discharge)
  0.5 * sum(log(w)) - 0.5 * n * log(mean(w * (discharge - mu)^2))
}
