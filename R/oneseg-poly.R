# Polynomial rating curve.

#' Fit polynomial rating curve
#'
#' @param q A vector of streamflow data.
#' @param h A vector of stage data.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param degree Polynomial degree.
#' @param wts_code Weighting scheme: `"none"`, `"spec"`, or `"prop"`.
#' @param wts Optional vector of weights when `wts_code = "spec"`.
#' @param wts_tol Convergence tolerance for proportional-weight iteration.
#' @param wts_maxiter Maximum iterations for proportional-weight fitting.
#' @return An rc_poly object.
#' @examples
#' fit <- rc_poly(q, h, data = thompson, degree = 2)
#' predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
#' @export
rc_poly <- function(
  q,
  h,
  data = NULL,
  degree = 2,
  wts_code = c("none", "spec", "prop"),
  wts = NULL,
  wts_tol = 1e-6,
  wts_maxiter = 100
) {
  # error checks and warnings
  # q and h may name columns of `data`, or be vectors
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  wts_code <- rlang::arg_match(wts_code)

  # remove missing observations
  qh <- tidyr::drop_na(data.frame(qobs = q, hobs = h))
  q <- qh$qobs
  h <- qh$hobs
  qh_fit <- data.frame(q = q, h = h)

  # lm model formula
  lm_modform <- "q ~ h"
  term <- "h"
  for (i in 2:degree) {
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
    if (!is.numeric(wts)) {
      stop("invalid values of specified weights")
    }
    mod_poly <- stats::lm(
      formula = lm_modform,
      weights = wts,
      data = cbind(qh_fit, wts = wts)
    )
    mod_sum <- summary(mod_poly)
  } else {
    # fit using proportional weights - start using lm with no weights
    mod_poly <- stats::lm(formula = lm_modform, data = qh_fit)
    qp <- stats::predict(mod_poly)
    wts <- 1 / qp^2
    coefs_old <- as.numeric(mod_poly$coefficients)

    # create model formula for nls fit
    nls_modform <- "q ~ b0 + b1*h"
    for (i in 2:degree) {
      nls_modform <- paste0(nls_modform, " + b", i, "*h^", i)
    }

    # create list of starting values
    startlist <- as.list(coefs_old)
    names(startlist) <- paste0("b", 0:degree)

    # iterate until coefficient values converge
    for (i in 1:wts_maxiter) {
      mod_poly <- stats::nls(
        formula = nls_modform,
        weights = wts,
        data = cbind(qh_fit, wts = wts),
        start = startlist
      )
      coefs <- stats::coefficients(mod_poly)
      max_change <- max(abs((coefs - coefs_old) / coefs_old))
      if (max_change < wts_tol) {
        break
      }
      coefs_old <- coefs
      qp <- stats::predict(mod_poly)
      wts <- 1 / qp^2
    }
    mod_sum <- summary(mod_poly)
  }

  if (wts_code == "prop") {
    mod_form <- nls_modform
  } else {
    mod_form <- lm_modform
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  outlist <- list(
    qh_obs = qh,
    formula = mod_form,
    wts_code = wts_code,
    weights = wts,
    pars = stats::coefficients(mod_sum),
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
  if (predlim && object$wts_code == "spec") {
    message("Note: prediction limits cannot be computed for specified weights")
  }
  if (is.null(hpred)) {
    hpred <- rc_hpred_grid(object)
  }
  checkmate::assert_numeric(hpred, min.len = 1, finite = TRUE)
  hpred_df <- data.frame(h = hpred)
  mod <- object[["model"]]
  wts_code <- object$wts_code
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
