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
    new = data.frame(h = hpred),
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
  qh <- object[["qh_obs"]]
  # the two-segment fits name the column "h", the one-segment fits "hobs"
  h <- if ("h" %in% names(qh)) qh[["h"]] else qh[["hobs"]]
  checkmate::assert_numeric(h, min.len = 1L, any.missing = FALSE)
  seq(min(h), max(h), length.out = n)
}
