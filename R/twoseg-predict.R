# Prediction, and delta-method limits.

#' Predict discharge with confidence / prediction limits from a two-segment fit
#'
#' Evaluates the fitted rating curve on a grid of stage values and, optionally,
#' attaches confidence and/or prediction limits. This is a thin dispatcher: it
#' validates the arguments common to every method, then hands off to one of the
#' three `*_limits_2seg()` functions, which all take the same
#' `(object, hpred, conflev, predlev)` arguments and all return the same
#' columns.
#'
#' @section Choosing a method:
#' The default, `"delta"`, is fast and is the classical choice, but **it is not
#' reliable near the breakpoint**. The mean function is not differentiable at
#' `h = k`, so the linearisation switches form there and the interval jumps: on
#' one fitted curve the band widened from 29 to 242 m^3 s^-1 across the
#' breakpoint. In a simulation study a nominal 95% delta interval covered the
#' true curve only about two-thirds of the time just above the breakpoint,
#' against roughly 97% for `"boot"`. Away from the breakpoint it behaves
#' normally.
#'
#' So: `"delta"` for a quick look or where the breakpoint is not of interest,
#' and `"boot"` where the interval matters.
#'
#' @param object An `rc_nls_2seg` fit (from [rc_nls_2seg()]).
#' @param ... Passed on to the chosen limits function.
#' @param hpred Stage values at which to return limits. Defaults to
#'   [rc_hpred_grid()]: 1000 points spanning the observed stage range.
#' @param conflev Confidence level for the mean-curve (confidence) interval, or
#'   `NULL` to omit it.
#' @param predlev Confidence level for the prediction interval, or `NULL` to
#'   omit it.
#' @param method Interval method, matched by [rlang::arg_match()]:
#'   \itemize{
#'     \item `"delta"` (the default): [delta_limits_2seg()], the linearised
#'       delta method. Fast, but fixes the breakpoint at `k-hat` and so jumps
#'       there; see the section above.
#'     \item `"boot"`: [boot_limits_2seg()], which resamples the gaugings and
#'       refits, and so does not rest on the asymptotic normal at all. Much the
#'       slowest, since it refits the model `B` times, and the most trustworthy
#'       at the breakpoint.
#'     \item `"sim"`: [sim_limits_2seg()], which draws parameters from their
#'       asymptotic normal distribution and pushes each through the model.
#'       Sensitive to the parameter scale; see `space` in [sim_limits_2seg()].
#'   }
#' @return A data frame (tibble if \pkg{tibble} is available) with column `h`,
#'   the fitted discharge `fit`, and, when requested, `ci_lwr`/`ci_upr` and
#'   `pi_lwr`/`pi_upr`. Which columns are present depends only on which of
#'   `conflev` and `predlev` were given -- never on the method or the
#'   weighting. Where a quantity cannot be computed (prediction limits under
#'   `"spec"` weights) the column is returned as `NA`.
#'
#'   Named for the object's class, `rc_nls_2seg`, so `predict(object)`
#'   dispatches here. The one-segment models define their own `predict.rc_nls`
#'   for the `rc_nls` class; the two do not collide.
#' @seealso [delta_limits_2seg()], [boot_limits_2seg()], [sim_limits_2seg()].
#' @examples
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   fit <- rc_nls_2seg(Q, H, data = sauze, kstart = 1)
#'   hp <- c(1, 1.5, 2, 4)
#'
#'   # the default, and fast
#'   predict(fit, hpred = hp, conflev = 0.95)
#'
#'   # slower, but does not assume the breakpoint is known. Compare the two
#'   # either side of the breakpoint, at about 1.85 m: the delta band jumps
#'   # there, the bootstrap band does not.
#'   predict(fit, hpred = hp, conflev = 0.95, method = "boot", B = 50)
#' }
#'
#' # the columns returned never depend on the model or the method, so results
#' # from different approaches stack directly
#' one <- rc_nls(q, h, data = thompson)
#' two <- rc_nls_2seg(q, h, data = thompson, wts_code = "prop", kstart = 2)
#' rbind(
#'   predict(one, hpred = 3, conflev = 0.95),
#'   predict(two, hpred = 3, conflev = 0.95)
#' )
#' @export
predict.rc_nls_2seg <- function(
  object,
  ...,
  hpred = NULL,
  conflev = NULL,
  predlev = NULL,
  method = c("delta", "boot", "sim")
) {
  method <- rlang::arg_match(method)
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  # hpred is deliberately NOT resolved here: every *_limits_2seg() function
  # defaults it the same way, so it resolves once, in whichever one runs.
  limits_fun <- switch(
    method,
    delta = delta_limits_2seg,
    boot = boot_limits_2seg,
    sim = sim_limits_2seg
  )
  limits_fun(
    object,
    hpred = hpred,
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
#'   \item `"none"` / `"spec"`: [investr::predFit()].
#'   \item `"prop"`: [nlspw_limits()], which accounts for the proportional
#'     error structure.
#' }
#' The linearisation holds the breakpoint fixed at `k-hat`, so the gradient
#' switches abruptly from one segment's to the other's as `h` crosses it and
#' the band jumps. [boot_limits_2seg()] makes no smoothness assumption and is
#' continuous there; prefer it when the interval near the breakpoint matters.
#'
#' Prediction limits are not identified under `"spec"` weights -- the
#' per-observation error variances are supplied, not estimated, so there is no
#' single scatter to add to the mean curve. Those columns come back `NA`.
#'
#' @param object An `rc_nls_2seg` fit (from [rc_nls_2seg()]).
#' @param hpred Stage values at which to return limits. Defaults to
#'   [rc_hpred_grid()]: 1000 points spanning the observed stage range.
#' @param ... Passed on to [investr::predFit()] for the `"none"`/`"spec"`
#'   cases.
#' @param conflev,predlev Levels for the confidence and prediction intervals,
#'   or `NULL` to omit either.
#' @return A data frame (tibble if \pkg{tibble} is available); see
#'   [predict.rc_nls_2seg()] for the columns.
#' @examples
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   fit <- rc_nls_2seg(Q, H, data = sauze, kstart = 1)
#'   delta_limits_2seg(fit, hpred = c(1, 2, 4), conflev = 0.95)
#' }
#' @export
delta_limits_2seg <- function(
  object,
  hpred = NULL,
  ...,
  conflev = NULL,
  predlev = NULL
) {
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1L, finite = TRUE)
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  if (predlim && object$wts_code == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  hpred_df <- data.frame(h = hpred)
  mod <- object[["model"]]
  wts_code <- object$wts_code
  # point predictions (fitted mean discharge) at the requested stages
  # yvec <- unname(stats::predict(mod, newdata = hpred_df, ...))
  yvec <- unname(stats::predict(mod, newdata = hpred_df))
  out_df <- data.frame(h = hpred, fit = yvec)
  # Confidence limits (uncertainty in the mean curve). Route by weighting:
  # unweighted/user-weighted -> delta method; proportional -> nlspw_limits().
  if (conflim) {
    if (wts_code == "none" || wts_code == "spec") {
      ci_mat <- investr::predFit(
        mod,
        newdata = hpred_df,
        interval = "confidence",
        level = conflev,
        ...
      )
    } else {
      ci_mat <- nlspw_limits(
        mod,
        type = "confidence",
        level = conflev,
        hpred = hpred
      )
    }
    ci_mat <- ci_mat[, c("lwr", "upr"), drop = FALSE]
    colnames(ci_mat) <- paste0("ci_", colnames(ci_mat))
    out_df <- cbind(out_df, as.data.frame(ci_mat))
  }
  # Prediction limits (mean uncertainty + observation scatter). Under "spec"
  # weights the observation variances are not recoverable, so the columns are
  # returned as NA rather than dropped -- the caller gets the same shape of
  # answer whatever the weighting, as the other three methods also do.
  if (predlim) {
    if (wts_code == "spec") {
      out_df$pi_lwr <- NA_real_
      out_df$pi_upr <- NA_real_
    } else {
      if (wts_code == "none") {
        pi_mat <- investr::predFit(
          mod,
          newdata = hpred_df,
          interval = "prediction",
          level = predlev,
          ...
        )
      } else {
        pi_mat <- nlspw_limits(
          mod,
          type = "prediction",
          level = predlev,
          hpred = hpred
        )
      }
      pi_mat <- pi_mat[, c("lwr", "upr"), drop = FALSE]
      colnames(pi_mat) <- paste0("pi_", colnames(pi_mat))
      out_df <- cbind(out_df, as.data.frame(pi_mat))
    }
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    out_df <- tibble::as_tibble(out_df)
  }
  out_df
}
