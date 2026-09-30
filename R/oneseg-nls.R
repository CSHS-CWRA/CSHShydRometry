# Power law fitted by NLS on the natural scale.

#' Fit rating curve using nls on untransformed data
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge`, `stage` and `wts`
#'   are evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param wts How the scatter of the gaugings is modelled: `wts_none()` (or
#'   `"none"`, the default), `wts_prop()` (or `"prop"`), or `wts_spec()`
#'   with the weights. See [wts].
#' @param control Settings for [stats::nls()], as from [stats::nls.control()].
#' @return An `rc_nls` object; see [rating_curve] for its contents.
#' @examples
#' fit <- rc_nls(discharge, stage, data = thompson)
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#'
#' # constant coefficient of variation instead of constant variance
#' fit_prop <- rc_nls(discharge, stage, data = thompson, wts = "prop")
#' predict(fit_prop, stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_nls <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  wts = wts_none(),
  control = stats::nls.control(maxiter = 1000, tol = 1e-6)
) {
  ## error checks and warnings
  # discharge, stage and wts may use columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  checkmate::assert_list(control, names = "named")
  # the weighting scheme; specified weights are evaluated in `data`, and
  # kept aligned with the gaugings that remain
  weighting <- rc_resolve_wts(
    wts,
    data,
    keep = stats::complete.cases(discharge, stage)
  )
  wts_code <- weighting$type
  wts <- weighting$values

  ## remove missing observations
  qh <- rc_complete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage

  ## generate starting values using lm on log-transformed data
  cstart <- min(stage) - 0.1 * (max(stage) - min(stage))
  start_lm <- stats::lm(log(discharge) ~ log(stage - cstart))
  astart <- unname(exp(start_lm$coefficients[1]))
  bstart <- unname(start_lm$coefficients[2])
  irls <- NULL

  # use nls to determine optimal parameters - no weights or specified weights
  if (wts_code != "prop") {
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
    # proportional weights, by iterative reweighting from the log-log curve
    fit_fun <- function(wts, start) {
      stats::nls(
        discharge ~ a * (stage - c)^b,
        data = data.frame(discharge = discharge, stage = stage, wts = wts),
        weights = wts,
        start = start,
        control = control
      )
    }
    res <- rc_irls(
      fit_fun,
      yp = astart * (stage - cstart)^bstart,
      start = list(a = astart, b = bstart, c = cstart),
      tol = weighting$tol,
      maxiter = weighting$maxiter
    )
    mod_nls <- res$model
    wts <- res$weights
    irls <- res$irls
  }
  mod_sum <- summary(mod_nls)
  coefs <- stats::coef(mod_nls)

  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    gaugings = qh,
    pars = list(a = coefs[["a"]], b = coefs[["b"]], c = coefs[["c"]]),
    settings = list(wts = weighting, control = control),
    weights = wts,
    irls = irls,
    rse = mod_sum$sigma,
    model = mod_nls
  )
  structure(outlist, class = c("rc_nls", "rating_curve"))
}


#' Predict method for rc_nls objects
#'
#' @inheritParams predict.rc_log_ols
#' @param object An rc_nls object.
#' @export
predict.rc_nls <- function(
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
    stage <- rc_stage_grid(object)
  }
  checkmate::assert_numeric(stage, min.len = 1, finite = TRUE)
  stage_df <- data.frame(stage = stage)
  mod <- object[["model"]]
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
  if (requireNamespace("tibble", quietly = TRUE)) {
    out_df <- tibble::as_tibble(out_df)
  }
  out_df
}
