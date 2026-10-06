# Prediction, and delta-method limits.

#' Predict discharge with confidence / prediction limits from a two-segment fit
#'
#' Evaluates the fitted rating curve at the stages given and, optionally,
#' attaches confidence and/or prediction limits, by the delta method or a
#' bootstrap.
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
#' @section The bootstrap:
#' `method = "boot"` resamples the gaugings with replacement, refits the
#' two-segment curve to each resample with the arguments the fit was made
#' with (repeating the search over starting breakpoints), and summarises the
#' spread of the refitted curves. Under [var_spec()] the variances are
#' resampled along with the gaugings. A resample fails if its fit errors or,
#' under [var_prop()], if its reweighting does not converge; failed resamples
#' are redrawn, up to `B * max_tries_factor` attempts in all.
#'
#' Confidence limits are the pointwise percentiles of the refitted curves.
#' Prediction limits add the scatter of a new gauging to the spread of the
#' refitted curves, using the t distribution on the fit's residual degrees of
#' freedom. Under [var_spec()] they are `NA`, as for the delta method. The
#' number of resamples that fitted is recorded in `attr(, "B_success")`.
#'
#' @section What the limits assume:
#' Prediction limits, from either method, assume the scatter of the
#' gaugings about the curve is normal, with the spread the weighting scheme
#' describes: if it is skewed or has heavy tails, they can miss, especially
#' at high levels such as 0.99. Delta-method confidence limits assume the
#' estimates are close to normal, which is approximate, and poor near the
#' breakpoint; bootstrap confidence limits do not.
#'
#' @param object An `rc_2seg_powerlaw` fit (from [rc_2seg_powerlaw()]).
#' @param ... For `method = "boot"`, its settings:
#'   \describe{
#'     \item{`B`}{The number of resamples that must fit successfully
#'       (default 1000).}
#'     \item{`seed`}{An optional random seed, for reproducible limits.}
#'     \item{`max_tries_factor`}{A cap on the attempts, `B *
#'       max_tries_factor` (default 3), so that the bootstrap ends even if
#'       many resamples fail.}
#'   }
#'   For `method = "delta"`, must be empty. Any other argument is an error.
#' @param new_stage Stages at which to return limits. Defaults to
#'   1000 points spanning the observed stage range.
#' @param conflev Confidence level for the mean-curve (confidence) interval, or
#'   `NULL` to omit it.
#' @param predlev Confidence level for the prediction interval, or `NULL` to
#'   omit it.
#' @param method Interval method:
#'   \itemize{
#'     \item `"delta"` (the default): the linearised delta method. Fast, but
#'       fixes the breakpoint at its estimate and so jumps there; see
#'       "Choosing a method".
#'     \item `"boot"`: a bootstrap, which resamples the gaugings and refits;
#'       see "The bootstrap". Much the slowest, since it refits the model `B`
#'       times, and the most trustworthy at the breakpoint.
#'   }
#' @return A tibble with column `stage`,
#'   the fitted discharge `fit`, and, when requested, `ci_lwr`/`ci_upr` and
#'   `pi_lwr`/`pi_upr`. Which columns are present depends only on which of
#'   `conflev` and `predlev` were given -- never on the method or the
#'   weighting. Where a quantity cannot be computed (prediction limits under
#'   `"spec"` weights) the column is returned as `NA`.
#'
#'   Named for the object's class, `rc_2seg_powerlaw`, so `predict(object)`
#'   dispatches here. The one-segment models define their own `predict.rc_powerlaw`
#'   for the `rc_powerlaw` class; the two do not collide.
#' @examples
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   fit <- rc_2seg_powerlaw(Q, H, data = sauze, kstart = 1)
#'   hp <- c(1, 1.5, 2, 4)
#'
#'   # the default, and fast
#'   predict(fit, new_stage = hp, conflev = 0.95)
#'
#'   # slower, but does not assume the breakpoint is known. Compare the two
#'   # either side of the breakpoint, at about 1.85 m: the delta band jumps
#'   # there, the bootstrap band does not.
#'   predict(fit, new_stage = hp, conflev = 0.95, method = "boot", B = 50)
#' }
#'
#' # the columns returned never depend on the model or the method, so results
#' # from different approaches stack directly
#' one <- rc_powerlaw(discharge, stage, data = thompson)
#' two <- rc_2seg_powerlaw(discharge, stage, data = thompson, variance = "prop", kstart = 2)
#' rbind(
#'   predict(one, new_stage = 3, conflev = 0.95),
#'   predict(two, new_stage = 3, conflev = 0.95)
#' )
#' @export
predict.rc_2seg_powerlaw <- function(
  object,
  ...,
  new_stage = NULL,
  conflev = NULL,
  predlev = NULL,
  method = c("delta", "boot")
) {
  method <- rlang::arg_match(method)
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  # the work is done by an internal function for each method, which resolves
  # new_stage and checks that ... holds only arguments it takes
  limits_fun <- switch(
    method,
    delta = delta_limits_2seg,
    boot = boot_limits_2seg
  )
  limits_fun(
    object,
    new_stage = new_stage,
    conflev = conflev,
    predlev = predlev,
    ...
  )
}


#' Delta-method limits for a two-segment fit
#'
#' The `method = "delta"` back end of [predict.rc_2seg_powerlaw()]. It
#' linearises the curve about the estimates, holding the breakpoint fixed, so
#' the band jumps at the breakpoint. Under `var_none()` and `var_spec()` it
#' uses [investr::predFit()]; under `var_prop()`, `nlspw_limits()`.
#'
#' @param object An `rc_2seg_powerlaw` fit.
#' @param ... Must be empty.
#' @param new_stage,conflev,predlev As in [predict.rc_2seg_powerlaw()].
#' @return A tibble; see [predict.rc_2seg_powerlaw()] for the columns.
#' @noRd
delta_limits_2seg <- function(
  object,
  ...,
  new_stage = NULL,
  conflev = NULL,
  predlev = NULL
) {
  rlang::check_dots_empty()
  stage <- if (is.null(new_stage)) stage_grid(object) else new_stage
  checkmate::assert_numeric(stage, min.len = 1L, finite = TRUE)
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  if (predlim && object$settings$variance$type == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
  wts_code <- object$settings$variance$type
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
        level = conflev
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
          level = predlev
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
