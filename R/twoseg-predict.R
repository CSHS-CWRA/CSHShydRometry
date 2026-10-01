# Prediction, and delta-method limits.

#' Predict discharge with confidence / prediction limits from a two-segment fit
#'
#' Evaluates the fitted rating curve on a grid of stage values and, optionally,
#' attaches confidence and/or prediction limits. This is a thin dispatcher: it
#' validates the arguments common to every method, then hands off to one of the
#' `*_limits_2seg()` functions, which all take the same
#' `(object, stage, conflev, predlev)` arguments and all return the same
#' columns.
#'
#' @section Choosing a method:
#' The default, `"delta"`, is fast and is the classical choice, but **it is not
#' reliable near the breakpoint**. The mean function is not differentiable at
#' `stage = k`, so the linearisation switches form there and the interval jumps: on
#' one fitted curve the band widened from 29 to 242 m^3 s^-1 across the
#' breakpoint. In a simulation study a nominal 95% delta interval covered the
#' true curve only about two-thirds of the time just above the breakpoint,
#' against roughly 97% for `"boot"`. Away from the breakpoint it behaves
#' normally.
#'
#' So: `"delta"` for a quick look or where the breakpoint is not of interest,
#' and `"boot"` where the interval matters.
#'
#' A third approach -- drawing parameter vectors from their asymptotic normal
#' distribution and pushing each draw through the model -- was tried and
#' dropped. When a segment is poorly identified, as the lower segment at Sauze
#' is, the draws too often land on impossible curves: negative or
#' astronomically large discharges, giving limits that are meaningless.
#'
#' @param object An `rc_2seg_nls` fit (from [rc_2seg_nls()]).
#' @param ... Passed on to the chosen limits function.
#' @param stage Stages at which to return limits. Defaults to
#'   1000 points spanning the observed stage range.
#' @param conflev Confidence level for the mean-curve (confidence) interval, or
#'   `NULL` to omit it.
#' @param predlev Confidence level for the prediction interval, or `NULL` to
#'   omit it.
#' @param method Interval method:
#'   \itemize{
#'     \item `"delta"` (the default): [delta_limits_2seg()], the linearised
#'       delta method. Fast, but fixes the breakpoint at `k-hat` and so jumps
#'       there; see the section above.
#'     \item `"boot"`: [boot_limits_2seg()], which resamples the gaugings and
#'       refits, and so does not rest on the asymptotic normal at all. Much the
#'       slowest, since it refits the model `B` times, and the most trustworthy
#'       at the breakpoint.
#'   }
#' @return A tibble with column `stage`,
#'   the fitted discharge `fit`, and, when requested, `ci_lwr`/`ci_upr` and
#'   `pi_lwr`/`pi_upr`. Which columns are present depends only on which of
#'   `conflev` and `predlev` were given -- never on the method or the
#'   weighting. Where a quantity cannot be computed (prediction limits under
#'   `"spec"` weights) the column is returned as `NA`.
#'
#'   Named for the object's class, `rc_2seg_nls`, so `predict(object)`
#'   dispatches here. The one-segment models define their own `predict.rc_nls`
#'   for the `rc_nls` class; the two do not collide.
#' @seealso [delta_limits_2seg()], [boot_limits_2seg()].
#' @examples
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   fit <- rc_2seg_nls(Q, H, data = sauze, kstart = 1)
#'   hp <- c(1, 1.5, 2, 4)
#'
#'   # the default, and fast
#'   predict(fit, stage = hp, conflev = 0.95)
#'
#'   # slower, but does not assume the breakpoint is known. Compare the two
#'   # either side of the breakpoint, at about 1.85 m: the delta band jumps
#'   # there, the bootstrap band does not.
#'   predict(fit, stage = hp, conflev = 0.95, method = "boot", B = 50)
#' }
#'
#' # the columns returned never depend on the model or the method, so results
#' # from different approaches stack directly
#' one <- rc_nls(discharge, stage, data = thompson)
#' two <- rc_2seg_nls(discharge, stage, data = thompson, wts = "prop", kstart = 2)
#' rbind(
#'   predict(one, stage = 3, conflev = 0.95),
#'   predict(two, stage = 3, conflev = 0.95)
#' )
#' @export
predict.rc_2seg_nls <- function(
  object,
  ...,
  stage = NULL,
  conflev = NULL,
  predlev = NULL,
  method = c("delta", "boot")
) {
  method <- rlang::arg_match(method)
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  # stage is deliberately NOT resolved here: every *_limits_2seg() function
  # defaults it the same way, so it resolves once, in whichever one runs.
  limits_fun <- switch(
    method,
    delta = delta_limits_2seg,
    boot = boot_limits_2seg
  )
  limits_fun(
    object,
    stage = stage,
    conflev = conflev,
    predlev = predlev,
    ...
  )
}


#' Delta-method confidence and prediction limits for a two-segment fit
#'
#' Linearises the fitted curve about `theta-hat` and propagates the parameter
#' covariance through that linearisation. Routed by weighting:
#' \itemize{
#'   \item [wts_none()] / [wts_spec()]: [investr::predFit()].
#'   \item [wts_prop()]: the delta method adapted to the proportional
#'     error structure.
#' }
#' The linearisation holds the breakpoint fixed at `k-hat`, so the gradient
#' switches abruptly from one segment's to the other's as `stage` crosses it and
#' the band jumps. [boot_limits_2seg()] makes no smoothness assumption and is
#' continuous there; prefer it when the interval near the breakpoint matters.
#'
#' Prediction limits are not identified under `"spec"` weights -- the
#' per-observation error variances are supplied, not estimated, so there is no
#' single scatter to add to the mean curve. Those columns come back `NA`.
#'
#' @param object An `rc_2seg_nls` fit (from [rc_2seg_nls()]).
#' @param stage Stages at which to return limits. Defaults to
#'   1000 points spanning the observed stage range.
#' @param ... Passed on to [investr::predFit()] for the [wts_none()]/[wts_spec()]
#'   cases.
#' @param conflev,predlev Levels for the confidence and prediction intervals,
#'   or `NULL` to omit either.
#' @return A tibble; see
#'   [predict.rc_2seg_nls()] for the columns.
#' @examples
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   fit <- rc_2seg_nls(Q, H, data = sauze, kstart = 1)
#'   delta_limits_2seg(fit, stage = c(1, 2, 4), conflev = 0.95)
#' }
#' @export
delta_limits_2seg <- function(
  object,
  ...,
  stage = NULL,
  conflev = NULL,
  predlev = NULL
) {
  if (is.null(stage)) {
    stage <- stage_grid(object)
  }
  checkmate::assert_numeric(stage, min.len = 1L, finite = TRUE)
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  if (predlim && object$settings$wts$type == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
  wts_code <- object$settings$wts$type
  # point predictions (fitted mean discharge) at the requested stages
  yvec <- unname(stats::predict(mod, newdata = stage_df))
  out_df <- data.frame(stage = stage, fit = yvec)
  # Confidence limits (uncertainty in the mean curve). Route by weighting:
  # unweighted/user-weighted -> delta method; proportional -> nlspw_limits().
  if (conflim) {
    if (wts_code == "none" || wts_code == "spec") {
      ci_mat <- investr::predFit(
        mod,
        newdata = stage_df,
        interval = "confidence",
        level = conflev,
        ...
      )
    } else {
      ci_mat <- nlspw_limits(
        mod,
        type = "confidence",
        level = conflev,
        stage = stage
      )
    }
    ci_mat <- ci_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(ci_mat) <- paste0("ci_", colnames(ci_mat))
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  # Prediction limits (mean uncertainty + observation scatter). Under "spec"
  # weights the observation variances are not recoverable, so the columns are
  # returned as NA rather than dropped -- the caller gets the same shape of
  # answer whatever the weighting, as the other methods also do.
  if (predlim) {
    if (wts_code == "spec") {
      out_df$pi_lwr <- NA_real_
      out_df$pi_upr <- NA_real_
    } else {
      if (wts_code == "none") {
        pi_mat <- investr::predFit(
          mod,
          newdata = stage_df,
          interval = "prediction",
          level = predlev,
          ...
        )
      } else {
        pi_mat <- nlspw_limits(
          mod,
          type = "prediction",
          level = predlev,
          stage = stage
        )
      }
      pi_mat <- pi_mat[, c("lwr", "upr"), drop = FALSE]
      colnames(pi_mat) <- paste0("pi_", colnames(pi_mat))
      out_df <- cbind(out_df, as.data.frame(pi_mat))
    }
  }
  out_df <- tibble::as_tibble(out_df)
  out_df
}
