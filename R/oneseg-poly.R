# Polynomial rating curve.

#' Fit polynomial rating curve
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge` and `stage` are
#'   evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param degree Polynomial degree; positive whole number. Default 2.
#' @param variance `r lifecycle::badge("experimental")` How the variance of
#'   the gaugings about the curve is modelled: `var_none()` (or `"none"`, the default), `var_prop()` (or
#'   `"prop"`), or `var_spec()` with the variances. See [variance].
#' @return An `rc_poly` object; see [rating_curve] for its contents. The
#'   coefficients are named `b0`, `b1`, ... by power of `stage`.
#' @examples
#' fit <- rc_poly(discharge, stage, data = thompson, degree = 2)
#' predict(fit, new_stage = c(1, 3, 6), conflev = 0.95)
#' @export
rc_poly <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  degree = 2,
  variance = var_none()
) {
  # error checks and warnings
  # discharge and stage may use columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  checkmate::assert_count(degree, positive = TRUE)
  # the weighting scheme; specified weights are kept aligned with the
  # gaugings that remain
  weighting <- resolve_variance(
    variance,
    keep = stats::complete.cases(discharge, stage),
    fitter = "rc_poly"
  )
  wts_code <- weighting$type
  # the engine fits with weights, the reciprocals of the variances
  wts <- if (weighting$type == "spec") 1 / weighting$values

  # remove missing observations
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  qh_fit <- data.frame(discharge = discharge, stage = stage)
  irls <- NULL

  # lm model formula
  lm_modform <- "discharge ~ stage"
  term <- "stage"
  for (i in seq_len(degree)[-1]) {
    term <- paste0(term, "*stage")
    lm_modform <- paste0(lm_modform, " + I(", term, ")")
  }

  # fit model
  if (wts_code == "none") {
    # use ols to determine optimal parameters
    mod_poly <- stats::lm(formula = lm_modform, data = qh_fit)
    mod_sum <- summary(mod_poly)
  } else if (wts_code == "spec") {
    # use specified weights
    checkmate::assert_numeric(wts, len = length(discharge), .var.name = "wts")
    mod_poly <- stats::lm(
      formula = lm_modform,
      weights = wts,
      data = cbind(qh_fit, wts = wts)
    )
    mod_sum <- summary(mod_poly)
  } else {
    # proportional weights, by iterative reweighting from the unweighted fit
    mod_ols <- stats::lm(formula = lm_modform, data = qh_fit)

    # create model formula for nls fit
    nls_modform <- "discharge ~ b0 + b1*stage"
    for (i in seq_len(degree)[-1]) {
      nls_modform <- paste0(nls_modform, " + b", i, "*stage^", i)
    }
    nls_modform <- stats::as.formula(nls_modform)

    # starting values from the unweighted fit
    startlist <- as.list(unname(stats::coef(mod_ols)))
    names(startlist) <- paste0("b", 0:degree)

    fit_fun <- function(wts, start) {
      stats::nls(
        formula = nls_modform,
        weights = wts,
        data = cbind(qh_fit, wts = wts),
        start = start
      )
    }
    res <- reweight_in_rounds(
      fit_fun,
      yp = as.numeric(stats::predict(mod_ols)),
      start = startlist,
      tol = weighting$tol,
      maxiter = weighting$maxiter
    )
    mod_poly <- res$model
    wts <- res$weights
    irls <- res$irls
    mod_sum <- summary(mod_poly)
  }
  if (wts_code == "none") {
    wts <- rep(1, length(discharge))
  }
  coefs <- unname(stats::coef(mod_poly))
  pars <- as.list(coefs)
  names(pars) <- paste0("b", 0:degree)
  qh <- tibble::as_tibble(qh)
  outlist <- list(
    gaugings = qh,
    curve_parameters = pars,
    settings = list(degree = degree, variance = weighting),
    weights_used = wts,
    irls = irls,
    variance = weighting,
    rse = mod_sum$sigma,
    model = mod_poly
  )
  structure(outlist, class = c("rc_poly", "rating_curve"))
}


#' Predict method for rc_poly objects
#'
#' @inheritParams predict.rc_powerlaw
#' @inheritSection predict.rc_powerlaw What the limits assume
#' @param object An rc_poly object.
#' @export
predict.rc_poly <- function(
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
  yvec <- unname(stats::predict(mod, newdata = stage_df))
  out_df <- data.frame(stage = stage, fit = yvec)
  if (conflim) {
    if (wts_code != "prop") {
      ci_mat <- stats::predict(
        mod,
        newdata = stage_df,
        interval = "confidence",
        level = conflev
      )
    } else {
      # proportional weights: the same limits as for rc_powerlaw()
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
  # depends only on what was asked for.
  if (predlim && wts_code == "spec") {
    out_df$pi_lwr <- NA_real_
    out_df$pi_upr <- NA_real_
  }
  if (predlim && wts_code != "spec") {
    if (wts_code == "none") {
      pi_mat <- stats::predict(
        mod,
        newdata = stage_df,
        interval = "prediction",
        level = predlev
      )
    } else if (wts_code == "prop") {
      # proportional weights: the same limits as for rc_powerlaw()
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
