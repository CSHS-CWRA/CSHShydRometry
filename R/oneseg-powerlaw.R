# Power law fitted by least squares on the original scale.

#' Fit a power-law rating curve on the original scale
#'
#' Fits the power law \eqn{Q = a (h - c)^b} by nonlinear least squares on the
#' original scale of discharge, with the scatter of the gaugings modelled as
#' set by `variance`. The curve estimates the mean discharge at each stage.
#'
#' @section The offset:
#' \eqn{c} is the offset: for a single power law, the stage at which the flow
#' would stop. It is estimated unless it is given, through `offset`, say from
#' a survey of the control. A given offset is held fixed, so the limits from
#' [predict()] carry no uncertainty in it. Fixing it at the value estimated by
#' a first fit treats an estimate as known: the curve is the same, but the
#' limits are narrower than they should be.
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge` and `stage` are
#'   evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param variance `r lifecycle::badge("experimental")` How the variance of
#'   the gaugings about the curve is modelled: `var_none()` (or `"none"`, the default), `var_prop()` (or
#'   `"prop"`), `var_power()` (or `"power"`), or `var_spec()` with the
#'   variances. See [variance].
#' @param offset The offset, \eqn{c} in the formula: `NULL` (the default) to
#'   estimate it, or a known value, below every gauged stage, to hold it
#'   fixed.
#' @param control Settings for [stats::nls()], as from [stats::nls.control()].
#' @return An `rc_powerlaw` object; see [rating_curve] for its contents.
#' @examples
#' fit <- rc_powerlaw(discharge, stage, data = thompson)
#' predict(fit, new_stage = c(1, 3, 6), conflev = 0.95)
#'
#' # constant coefficient of variation instead of constant variance
#' fit_prop <- rc_powerlaw(discharge, stage, data = thompson, variance = "prop")
#' predict(fit_prop, new_stage = c(1, 3, 6), conflev = 0.95)
#'
#' # the offset known, say from a survey of the control
#' rc_powerlaw(discharge, stage, data = thompson, variance = "prop", offset = -1.3)
#' @export
rc_powerlaw <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  variance = var_none(),
  offset = NULL,
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
  checkmate::assert_number(offset, null.ok = TRUE, finite = TRUE)
  checkmate::assert_list(control, names = "named")
  # the weighting scheme; specified weights are kept aligned with the
  # gaugings that remain
  weighting <- resolve_variance(
    variance,
    keep = stats::complete.cases(discharge, stage),
    fitter = "rc_powerlaw",
    power_ok = TRUE
  )
  # as given, for refitting; `weighting` may gain an estimate below
  wts_given <- weighting
  wts_code <- weighting$type
  # the engine fits with weights, the reciprocals of the variances
  wts <- if (weighting$type == "spec") 1 / weighting$values

  ## remove missing observations
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  if (!is.null(offset) && offset >= min(stage)) {
    stop("`offset` must be below every gauged stage.", call. = FALSE)
  }

  ## the curve; a given offset enters it as a number, not a parameter
  formula <- if (is.null(offset)) {
    discharge ~ a * (stage - c)^b
  } else {
    eval(bquote(discharge ~ a * (stage - .(offset))^b))
  }

  ## generate starting values using lm on log-transformed data
  cstart <- if (is.null(offset)) {
    min(stage) - 0.1 * (max(stage) - min(stage))
  } else {
    offset
  }
  start_lm <- stats::lm(log(discharge) ~ log(stage - cstart))
  astart <- unname(exp(start_lm$coefficients[1]))
  bstart <- unname(start_lm$coefficients[2])
  start <- list(a = astart, b = bstart)
  if (is.null(offset)) {
    start$c <- cstart
  }
  irls <- NULL

  # use nls to determine optimal parameters - no weights or specified weights
  if (wts_code %in% c("none", "spec")) {
    if (wts_code == "none") {
      wts <- rep(1, length(discharge))
    }
    checkmate::assert_numeric(wts, len = length(discharge), .var.name = "wts")
    mod_nls <- stats::nls(
      formula,
      data = data.frame(discharge = discharge, stage = stage),
      weights = wts,
      start = start,
      control = control
    )
  } else {
    # proportional weights, by iterative reweighting from the log-log curve;
    # under var_power() this is the starting point for estimating the power
    rounds <- if (wts_code == "prop") weighting else var_prop()
    fit_fun <- function(wts, start) {
      stats::nls(
        formula,
        data = data.frame(discharge = discharge, stage = stage, wts = wts),
        weights = wts,
        start = start,
        control = control
      )
    }
    res <- reweight_in_rounds(
      fit_fun,
      yp = astart * (stage - cstart)^bstart,
      start = start,
      tol = rounds$tol,
      maxiter = rounds$maxiter
    )
    mod_nls <- res$model
    wts <- res$weights
    irls <- res$irls
    if (wts_code == "power") {
      # estimate the power with the curve, by generalised least squares
      mod_nls <- nlme::gnls(
        formula,
        data = data.frame(discharge = discharge, stage = stage),
        start = as.list(stats::coef(mod_nls)),
        weights = weighting$varfunc,
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
  c <- if (is.null(offset)) coefs[["c"]] else offset

  qh <- tibble::as_tibble(qh)
  outlist <- list(
    gaugings = qh,
    pars = list(a = coefs[["a"]], b = coefs[["b"]], c = c),
    settings = list(variance = wts_given, offset = offset, control = control),
    weights_used = wts,
    irls = irls,
    variance = weighting,
    rse = mod_sum$sigma,
    model = mod_nls
  )
  structure(outlist, class = c("rc_powerlaw", "rating_curve"))
}


#' Predict method for rc_powerlaw objects
#'
#' @param object An `rc_powerlaw` object.
#' @param new_stage Stages at which to predict discharge. If `NULL`, a grid of
#'   1000 equally spaced stages spanning the observed range is used.
#' @param ... Must be empty. Present so that every argument after it has to
#'   be named in full; a misspelt or unsupported argument is an error, not
#'   silently ignored.
#' @param conflev The confidence level for the confidence limits; if NULL, no
#'   confidence limits are returned.
#' @param predlev The prediction level for the prediction limits; if NULL, no
#'   prediction limits are returned.
#' @section What the limits assume:
#' The limits use the t distribution, so they assume the scatter of the
#' gaugings about the curve is normal (on the log scale, for
#' [rc_powerlaw_log()]), with the spread the weighting scheme describes.
#' Confidence limits depend on this only mildly, because estimates average
#' over many gaugings. Prediction limits depend on it directly: if the
#' scatter is skewed or has heavy tails, they can miss, especially at high
#' levels such as 0.99. Check the residuals before relying on them.
#'
#' Even when the scatter is normal, most limits are approximations: a 95%
#' interval covers the truth roughly, not exactly, 95% of the time. That is
#' because the curve is nonlinear in its parameters and is approximated by
#' a straight line about the estimates (the delta method), because the
#' weights are themselves estimated from the fit (under [var_prop()] and
#' [var_power()]), or, for [rc_loess()], because of the smoother's
#' approximations. The approximation is good with plenty of gaugings, and
#' poorer with few, or beyond the range of the gaugings. The limits are
#' exact only for fits that are linear in their parameters with weights
#' fixed in advance: [rc_poly()] with [var_none()] or [var_spec()], and
#' [rc_powerlaw_log()] with `offset` given.
#' @return A tibble with the stages (`stage`), the predicted discharges
#'   (`fit`) and, if requested, the lower and upper confidence limits
#'   (`ci_lwr`, `ci_upr`) and prediction limits (`pi_lwr`, `pi_upr`).
#' @export
predict.rc_powerlaw <- function(
  object,
  ...,
  new_stage = NULL,
  conflev = NULL,
  predlev = NULL
) {
  rlang::check_dots_empty()
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  wts_code <- object$settings$variance$type
  if (predlim && wts_code == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  stage <- if (is.null(new_stage)) stage_grid(object) else new_stage
  checkmate::assert_numeric(stage, min.len = 1, finite = TRUE)
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
  if (inherits(mod, "gnls")) {
    # the power was estimated, by nlme::gnls()
    return(gnls_limits(
      mod, stage_df, object$pars, object$variance$exponent, conflev, predlev
    ))
  }
  yvec <- unname(stats::predict(mod, newdata = stage_df))
  out_df <- data.frame(stage = stage, fit = yvec)
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
  out_df <- tibble::as_tibble(out_df)
  out_df
}


#' Limits for a power law fitted with the power of its scatter estimated
#'
#' Delta-method limits for a power law fitted by [nlme::gnls()] with the
#' scatter proportional to a power of the fitted flow. The curve's standard
#' error is `sqrt(g' V g)`, where `g` is the gradient of `a (h - c)^b` in its
#' estimated parameters (`c` may be fixed) and `V` their covariance matrix. A new gauging adds its scatter,
#' `sigma^2 fit^(2 power)`. The estimated power is treated as known, and the
#' limits use the t distribution on the fit's residual degrees of freedom.
#'
#' @param mod An [nlme::gnls()] fit of `a * (stage - c)^b`.
#' @param stage_df Data frame of stages, with column `stage`.
#' @param pars The curve's parameters, `a`, `b` and `c`.
#' @param power The estimated power of the scatter.
#' @param conflev,predlev Levels, or `NULL` to omit those limits.
#' @return A tibble of predictions and limits.
#' @noRd
gnls_limits <- function(mod, stage_df, pars, power, conflev, predlev) {
  a <- pars[["a"]]
  b <- pars[["b"]]
  c <- pars[["c"]]
  depth <- stage_df$stage - c
  fit <- a * depth^b
  gradient <- cbind(
    a = depth^b,
    b = fit * log(depth),
    c = -a * b * depth^(b - 1)
  )
  estimated <- names(stats::coef(mod))
  gradient <- gradient[, estimated, drop = FALSE]
  vcov_mat <- stats::vcov(mod)[estimated, estimated, drop = FALSE]
  se_fit <- sqrt(rowSums((gradient %*% vcov_mat) * gradient))
  df <- mod$dims$N - mod$dims$p
  out_df <- data.frame(stage = stage_df$stage, fit = fit)
  if (!is.null(conflev)) {
    margin <- stats::qt((1 + conflev) / 2, df) * se_fit
    out_df$ci_lwr <- fit - margin
    out_df$ci_upr <- fit + margin
  }
  if (!is.null(predlev)) {
    se_new <- sqrt(se_fit^2 + mod$sigma^2 * fit^(2 * power))
    margin <- stats::qt((1 + predlev) / 2, df) * se_new
    out_df$pi_lwr <- fit - margin
    out_df$pi_upr <- fit + margin
  }
  tibble::as_tibble(out_df)
}
