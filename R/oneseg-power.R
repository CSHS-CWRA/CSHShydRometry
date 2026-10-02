# Power law fitted by least squares on the original scale.

#' Fit a power-law rating curve on the original scale
#'
#' Fits the power law \eqn{Q = a (h - c)^b} by nonlinear least squares on the
#' original scale of discharge, with the scatter of the gaugings modelled as
#' set by `wts`. The curve estimates the mean discharge at each stage.
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge` and `stage` are
#'   evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param wts How the scatter of the gaugings is modelled: `wts_none()` (or
#'   `"none"`, the default), `wts_prop()` (or `"prop"`), or `wts_spec()`
#'   with the weights. See [wts].
#' @param control Settings for [stats::nls()], as from [stats::nls.control()].
#' @return An `rc_power` object; see [rating_curve] for its contents.
#' @examples
#' fit <- rc_power(discharge, stage, data = thompson)
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#'
#' # constant coefficient of variation instead of constant variance
#' fit_prop <- rc_power(discharge, stage, data = thompson, wts = "prop")
#' predict(fit_prop, stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_power <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  wts = wts_none(),
  control = stats::nls.control(maxiter = 1000, tol = 1e-6)
) {
  ## error checks and warnings
  # discharge and stage may use columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  checkmate::assert_list(control, names = "named")
  # the weighting scheme; specified weights are kept aligned with the
  # gaugings that remain
  weighting <- resolve_wts(
    wts,
    keep = stats::complete.cases(discharge, stage),
    fitter = "rc_power",
    power_ok = TRUE
  )
  # as given, for refitting; `weighting` may gain an estimate below
  wts_given <- weighting
  wts_code <- weighting$type
  wts <- weighting$values

  ## remove missing observations
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage

  ## generate starting values using lm on log-transformed data
  cstart <- min(stage) - 0.1 * (max(stage) - min(stage))
  start_lm <- stats::lm(log(discharge) ~ log(stage - cstart))
  astart <- unname(exp(start_lm$coefficients[1]))
  bstart <- unname(start_lm$coefficients[2])
  irls <- NULL

  # use nls to determine optimal parameters - no weights or specified weights
  if (wts_code %in% c("none", "spec")) {
    if (wts_code == "none") {
      wts <- rep(1, length(discharge))
    }
    checkmate::assert_numeric(wts, len = length(discharge), .var.name = "wts")
    mod_nls <- stats::nls(
      discharge ~ a * (stage - c)^b,
      data = data.frame(discharge = discharge, stage = stage),
      weights = wts,
      start = list(a = astart, b = bstart, c = cstart),
      control = control
    )
  } else {
    # proportional weights, by iterative reweighting from the log-log curve;
    # under wts_power() this is the starting point for estimating the power
    rounds <- if (wts_code == "prop") weighting else wts_prop()
    fit_fun <- function(wts, start) {
      stats::nls(
        discharge ~ a * (stage - c)^b,
        data = data.frame(discharge = discharge, stage = stage, wts = wts),
        weights = wts,
        start = start,
        control = control
      )
    }
    res <- reweight_in_rounds(
      fit_fun,
      yp = astart * (stage - cstart)^bstart,
      start = list(a = astart, b = bstart, c = cstart),
      tol = rounds$tol,
      maxiter = rounds$maxiter
    )
    mod_nls <- res$model
    wts <- res$weights
    irls <- res$irls
    if (wts_code == "power") {
      # estimate the power with the curve, by generalised least squares
      mod_nls <- nlme::gnls(
        discharge ~ a * (stage - c)^b,
        data = data.frame(discharge = discharge, stage = stage),
        start = as.list(stats::coef(mod_nls)),
        weights = weighting$variance,
        control = nlme::gnlsControl(maxIter = 1e5, minScale = 1e-5)
      )
      weighting$exponent <- unname(stats::coef(
        mod_nls$modelStruct$varStruct,
        unconstrained = FALSE
      ))
      wts <- 1 / as.numeric(stats::fitted(mod_nls))^(2 * weighting$exponent)
      irls <- NULL
    }
  }
  mod_sum <- summary(mod_nls)
  coefs <- stats::coef(mod_nls)

  qh <- tibble::as_tibble(qh)
  outlist <- list(
    gaugings = qh,
    pars = list(a = coefs[["a"]], b = coefs[["b"]], c = coefs[["c"]]),
    settings = list(wts = wts_given, control = control),
    weights_used = wts,
    irls = irls,
    wts = weighting,
    rse = mod_sum$sigma,
    model = mod_nls
  )
  structure(outlist, class = c("rc_power", "rating_curve"))
}


#' Predict method for rc_power objects
#'
#' @param object An `rc_power` object.
#' @param stage Stages at which to predict discharge. If `NULL`, a grid of
#'   1000 equally spaced stages spanning the observed range is used.
#' @param ... Additional arguments passed to the inner `stats::predict()` function.
#' @param conflev The confidence level for the confidence limits; if NULL, no
#'   confidence limits are returned.
#' @param predlev The prediction level for the prediction limits; if NULL, no
#'   prediction limits are returned.
#' @return A tibble with the stages (`stage`), the predicted discharges
#'   (`fit`) and, if requested, the lower and upper confidence limits
#'   (`ci_lwr`, `ci_upr`) and prediction limits (`pi_lwr`, `pi_upr`).
#' @export
predict.rc_power <- function(
  object,
  ...,
  stage = NULL,
  conflev = NULL,
  predlev = NULL
) {
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  wts_code <- object$settings$wts$type
  if (predlim && wts_code == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  if (is.null(stage)) {
    stage <- stage_grid(object)
  }
  checkmate::assert_numeric(stage, min.len = 1, finite = TRUE)
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
  if (inherits(mod, "gnls")) {
    # the power was estimated: nlraa handles the gnls fit
    return(gnls_limits(mod, stage_df, conflev, predlev, ...))
  }
  yvec <- unname(stats::predict(mod, newdata = stage_df, ...))
  out_df <- data.frame(stage = stage, fit = yvec)
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
  # Under "spec" weights the observation variances are supplied rather than
  # estimated, so a new observation's scatter is not identified. Return the
  # columns as NA rather than dropping them, so that the set of columns
  # depends only on what was asked for (Dan, 2026-08-03; Paul agreed).
  if (predlim && wts_code == "spec") {
    out_df$pi_lwr <- NA_real_
    out_df$pi_upr <- NA_real_
  }
  if (predlim && wts_code != "spec") {
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
  out_df <- tibble::as_tibble(out_df)
  out_df
}


#' Limits for a power law fitted with the power of its scatter estimated
#'
#' @param mod An [nlme::gnls()] fit.
#' @param stage_df Data frame of stages, with column `stage`.
#' @param conflev,predlev Levels, or `NULL` to omit those limits.
#' @param ... Passed on to [nlraa::predict_gnls()].
#' @return A tibble of predictions and limits.
#' @noRd
gnls_limits <- function(mod, stage_df, conflev, predlev, ...) {
  limits <- function(interval, level) {
    lims <- nlraa::predict_gnls(
      mod,
      newdata = stage_df,
      interval = interval,
      level = level,
      ...
    )
    # predict_gnls() names the limits after their quantiles, e.g. Q2.5 and
    # Q97.5 at level 0.95; build the names the same way, and select by them
    lwr <- paste0("Q", (1 - level) / 2 * 100)
    upr <- paste0("Q", (1 - (1 - level) / 2) * 100)
    if (!all(c(lwr, upr) %in% colnames(lims))) {
      stop(
        "Unexpected columns from nlraa::predict_gnls(): expected ", lwr,
        " and ", upr, ", found ", paste(colnames(lims), collapse = ", "), ".",
        call. = FALSE
      )
    }
    list(lwr = unname(lims[, lwr]), upr = unname(lims[, upr]))
  }
  out_df <- data.frame(
    stage = stage_df$stage,
    fit = unname(nlraa::predict_gnls(mod, newdata = stage_df, ...))
  )
  if (!is.null(conflev)) {
    ci <- limits("confidence", conflev)
    out_df$ci_lwr <- ci$lwr
    out_df$ci_upr <- ci$upr
  }
  if (!is.null(predlev)) {
    pi <- limits("prediction", predlev)
    out_df$pi_lwr <- pi$lwr
    out_df$pi_upr <- pi$upr
  }
  tibble::as_tibble(out_df)
}
