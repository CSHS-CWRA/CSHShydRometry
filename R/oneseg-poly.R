# Polynomial rating curve.

#' Fit polynomial rating curve
#'
#' @param q A vector of streamflow data.
#' @param h A vector of stage data.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named in full.
#' @param degree Polynomial degree.
#' @param wts_code Weighting scheme: `"none"`, `"spec"`, or `"prop"`.
#' @param wts Optional vector of weights when `wts_code = "spec"`.
#' @param wts_tol Convergence tolerance under `wts_code = "prop"`: the
#'   reweighting stops once no fitted discharge changes by more than this
#'   fraction from one round to the next.
#' @param wts_maxiter Maximum number of reweighting rounds under
#'   `wts_code = "prop"`. Reaching it gives a warning.
#' @return An `rc_poly` object; see [rating_curve] for its contents. The
#'   coefficients are named `b0`, `b1`, ... by power of `h`.
#' @examples
#' fit <- rc_poly(q, h, data = thompson, degree = 2)
#' predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
#' @export
rc_poly <- function(
  q,
  h,
  ...,
  data = NULL,
  degree = 2,
  wts_code = c("none", "spec", "prop"),
  wts = NULL,
  wts_tol = 1e-6,
  wts_maxiter = 100
) {
  # error checks and warnings
  # q and h may name columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  checkmate::assert_count(degree, positive = TRUE)
  wts_code <- rlang::arg_match(wts_code)

  # remove missing observations
  # keep user-supplied weights aligned with the gaugings that remain
  if (length(wts) == length(q)) {
    wts <- wts[stats::complete.cases(q, h)]
  }
  qh <- rc_complete(q, h)
  q <- qh$q
  h <- qh$h
  qh_fit <- data.frame(q = q, h = h)
  wts_input <- wts
  irls <- NULL

  # lm model formula
  lm_modform <- "q ~ h"
  term <- "h"
  for (i in seq_len(degree)[-1]) {
    term <- paste0(term, "*h")
    lm_modform <- paste0(lm_modform, " + I(", term, ")")
  }

  # fit model
  if (wts_code == "none") {
    # use ols to determine optimal parameters
    mod_poly <- stats::lm(formula = lm_modform, data = qh_fit)
    mod_sum <- summary(mod_poly)
  } else if (wts_code == "spec") {
    # use specified weights
    checkmate::assert_numeric(wts, len = length(q), .var.name = "wts")
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
    nls_modform <- "q ~ b0 + b1*h"
    for (i in seq_len(degree)[-1]) {
      nls_modform <- paste0(nls_modform, " + b", i, "*h^", i)
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
    res <- rc_irls(
      fit_fun,
      yp = as.numeric(stats::predict(mod_ols)),
      start = startlist,
      wts_tol = wts_tol,
      wts_maxiter = wts_maxiter
    )
    mod_poly <- res$model
    wts <- res$weights
    irls <- res$irls
    mod_sum <- summary(mod_poly)
  }
  if (wts_code == "none") {
    wts <- rep(1, length(q))
  }
  coefs <- unname(stats::coef(mod_poly))
  pars <- as.list(coefs)
  names(pars) <- paste0("b", 0:degree)
  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    qh_obs = qh,
    pars = pars,
    settings = list(
      degree = degree,
      wts_code = wts_code,
      wts = wts_input,
      wts_tol = wts_tol,
      wts_maxiter = wts_maxiter
    ),
    weights = wts,
    irls = irls,
    rse = mod_sum$sigma,
    model = mod_poly
  )
  structure(outlist, class = c("rc_poly", "rating_curve"))
}


#' Predict method for rc_poly objects
#'
#' @inheritParams predict.rc_log_ols
#' @param object An rc_poly object.
#' @export
predict.rc_poly <- function(
  object,
  ...,
  hpred = NULL,
  conflev = NULL,
  predlev = NULL
) {
  checkmate::assert_number(conflev, null.ok = TRUE, lower = 0, upper = 1)
  checkmate::assert_number(predlev, null.ok = TRUE, lower = 0, upper = 1)
  predlim <- !is.null(predlev)
  conflim <- !is.null(conflev)
  wts_code <- object$settings$wts_code
  if (predlim && wts_code == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1, finite = TRUE)
  hpred_df <- data.frame(h = hpred)
  mod <- object[["model"]]
  yvec <- unname(stats::predict(mod, newdata = hpred_df, ...))
  out_df <- data.frame(h = hpred, fit = yvec)
  if (conflim) {
    if (wts_code != "prop") {
      ci_mat <- stats::predict(
        mod,
        newdata = hpred_df,
        interval = "confidence",
        level = conflev,
        ...
      )
    } else {
      # confidence limits, proportional weights
      cl_poly <- investr::predFit(
        mod,
        newdata = hpred_df,
        se.fit = TRUE
      )
      se_fit <- cl_poly$se.fit
      fit <- cl_poly$fit
      df <- cl_poly$df
      tc <- stats::qt(p = 0.5 + 0.5 * conflev, df = df)
      cl_lims <- data.frame(
        fit = fit,
        lwr = fit - tc * se_fit,
        upr = fit + tc * se_fit
      )
      ci_mat <- cl_lims
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
      pi_mat <- stats::predict(
        mod,
        newdata = hpred_df,
        interval = "prediction",
        level = predlev,
        ...
      )
    } else if (wts_code == "prop") {
      # compute weights for new observations if wts_code == "prop"
      qp <- stats::predict(mod, newdata = hpred_df)
      wtsp <- 1 / qp^2
      # prediction limits, proportional weights
      pl_poly <- investr::predFit(
        mod,
        newdata = hpred_df,
        se.fit = TRUE
      )
      se_fit <- pl_poly$se.fit
      res_scale <- pl_poly$residual.scale
      fit <- pl_poly$fit
      df <- pl_poly$df
      tc <- stats::qt(p = 0.5 + 0.5 * predlev, df = df)
      sp <- sqrt(se_fit^2 + (res_scale^2 / wtsp))
      pl_lims <- data.frame(
        fit = fit,
        lwr = fit - tc * sp,
        upr = fit + tc * sp
      )
      pi_mat <- pl_lims
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
