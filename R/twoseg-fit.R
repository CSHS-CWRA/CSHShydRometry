# Fitting the two-segment power-law rating curve.

#' Fit two-segment power-law rating curve using nls on untransformed data
#'
#' @param discharge <[`data-masking`][rlang::args_data_masking]> Discharge: a
#'   vector, or an expression evaluated in `data`, such as a column name.
#' @param stage <[`data-masking`][rlang::args_data_masking]> Stage: a vector,
#'   or an expression evaluated in `data`, such as a column name.
#' @param data Optional data frame in which `discharge`, `stage` and `wts`
#'   are evaluated.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named, which keeps calls readable and guards against
#'   positional mistakes.
#' @param controls How the two hydraulic controls combine above the
#'   breakpoint. `"successive"` (the default): the upper power law takes over
#'   from the lower one, with its coefficient set so that the two meet at the
#'   breakpoint. `"additive"`: the upper power law adds to the discharge the
#'   lower one carries at the breakpoint, as when flow spills onto a
#'   floodplain.
#' @param kstart Starting value(s) for the breakpoint `k`. `NULL`, the
#'   default, tries 10 values spread evenly across the search range; a
#'   numeric vector tries each of its values. See Details.
#' @param kfixed If `TRUE`, hold the breakpoint `k` fixed at `kstart`.
#' @param kbounds Lower and upper bounds for `k`, or `NULL` to keep at least
#'   three gaugings in each segment.
#' @param wts How the scatter of the gaugings is modelled: `wts_none()` (or
#'   `"none"`, the default), `wts_prop()` (or `"prop"`), or `wts_spec()`
#'   with the weights. See [wts].
#' @param control Settings for [stats::nls()], as from [stats::nls.control()].
#'   The fit uses the `"port"` algorithm, which ignores `tol`; its own
#'   settings, such as `rel.tol`, can be added to the list (see
#'   [stats::nls()]).
#'
#' @details
#' Fitting proceeds in three steps, the last two repeated for each starting
#' breakpoint:
#' \enumerate{
#'   \item Choose the breakpoint search range (from `kbounds`, or a default
#'         that leaves at least three gaugings in each segment) and the
#'         starting breakpoints to try.
#'   \item Derive starting values for the segment parameters by splitting the
#'         data at the starting breakpoint and fitting a linear model to each
#'         segment on the log-log scale (`log(discharge) ~ log(stage - c)`);
#'         `a = exp(intercept)`, `b = slope`.
#'   \item Fit all parameters jointly with [stats::nls()] using the "port"
#'         algorithm (which supports the box constraints in `lower`/`upper`).
#'         Under `wts_prop()` this fit is repeated, updating the weights
#'         from the current fitted values and starting from the previous
#'         estimates, until the fitted discharges stabilise.
#' }
#'
#' The two-segment fit is sensitive to where it starts: from some starting
#' breakpoints `nls()` fails outright, and from others it settles on a
#' local optimum. Trying several starts guards against both. Of the fits that
#' succeed, the one kept has the highest log-likelihood under the error model
#' (for `"none"` and `"spec"`, the smallest weighted residual sum of squares;
#' for `"prop"`, where the weights depend on the fit, the normal likelihood
#' with standard deviation proportional to the mean). Under `"prop"`, fits
#' whose reweighting converged are preferred. What each start led to is
#' recorded in `kstart_search`.
#'
#' Each start costs a full fit, and [boot_limits_2seg()] repeats the same
#' search for every resample. That is deliberate: a resample's best fit is
#' often found from a different start than the original's, so starting the
#' resamples only from the original estimate, or only from the optima the
#' original search found, would understate the uncertainty in the breakpoint.
#' A single `kstart` makes both the fit and the bootstrap faster, but both
#' then rest on that one start.
#'
#' @return An object of class `c("rc_2seg_power", "rating_curve")`; see
#'   [rating_curve] for its contents. `pars` holds the estimated parameters by
#'   type, one value per segment: under `"successive"`, `a` has a single value
#'   because the upper segment's coefficient is fixed by continuity, and under
#'   `"additive"`, `c` has a single value because the upper segment is
#'   measured from `k`. `coef()` gives the same estimates by their model
#'   names (`a1`, `b1`, `c1`, ...). `kstart_search` is a tibble with a row
#'   per starting breakpoint tried: the estimated breakpoint `k` it led to and
#'   the log-likelihood `loglik` of that fit, both `NA` where the fit failed.
#' @examples
#' # The Thompson is close to a single control, so its two-segment fit needs
#' # proportional weights to converge.
#' fit <- rc_2seg_power(discharge, stage, data = thompson, wts = "prop")
#' fit
#' coef(fit)
#'
#' # what each starting breakpoint led to
#' fit$kstart_search
#'
#' predict(fit, stage = c(1, 3, 6), conflev = 0.95)
#'
#' # A river with a clearer change of control is far less fussy. The Ardeche
#' # at Sauze, in the RBaM package, fits however the controls combine, and
#' # carries a reported uncertainty for every gauging.
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   succ <- rc_2seg_power(Q, H, data = sauze, kstart = 1)
#'   add <- rc_2seg_power(Q, H, data = sauze, controls = "additive", kstart = 1)
#'   c(successive = succ$pars[["k"]], additive = add$pars[["k"]])
#'
#'   # weights from the reported gauging uncertainties
#'   rc_2seg_power(
#'     Q,
#'     H,
#'     data = sauze,
#'     wts = wts_spec(1 / uQ^2),
#'     kstart = 1
#'   )
#' }
#' @export
rc_2seg_power <- function(
  discharge,
  stage,
  ...,
  data = NULL,
  controls = c("successive", "additive"),
  kstart = NULL,
  kfixed = FALSE,
  kbounds = NULL,
  wts = wts_none(),
  control = stats::nls.control(maxiter = 1000)
) {
  # -- 1. Inputs: tidy evaluation, checks, missing values ----
  # discharge, stage and wts may use columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  discharge <- rlang::eval_tidy(rlang::enquo(discharge), data)
  stage <- rlang::eval_tidy(rlang::enquo(stage), data)

  # error checks
  controls <- rlang::arg_match(controls)
  checkmate::assert_numeric(discharge, min.len = 1L)
  checkmate::assert_numeric(stage, len = length(discharge))
  checkmate::assert_flag(kfixed)
  checkmate::assert_numeric(kstart, min.len = 1L, finite = TRUE, null.ok = TRUE)
  checkmate::assert_numeric(kbounds, len = 2L, null.ok = TRUE)
  checkmate::assert_list(control, names = "named")
  # the weighting scheme; specified weights are evaluated in `data`, and
  # kept aligned with the gaugings that remain
  weighting <- resolve_wts(
    wts,
    data,
    keep = stats::complete.cases(discharge, stage),
    fitter = "rc_2seg_power"
  )
  wts_code <- weighting$type
  wts <- weighting$values
  # remove missing values, and check the number of observations
  qh <- drop_incomplete(discharge, stage)
  discharge <- qh$discharge
  stage <- qh$stage
  hsort <- sort(stage)
  n <- length(stage)
  if (n < 7) {
    stop("n < 7 - too few data points to fit two-segment curve")
  }
  # -- 2. Breakpoint: search range and starting values ----
  # The search range (klwr, kupr) comes from kbounds, or by default keeps at
  # least 3 observations in each segment (hsort[3] .. hsort[n-2]). The
  # starting values are those supplied, or a grid across the search range.
  kstart_input <- kstart
  if (kfixed) {
    # hold k fixed at the supplied value (lower == upper bound)
    if (length(kstart) != 1L) {
      stop("`kstart` must be a single value when `kfixed = TRUE`")
    }
    klwr <- kstart
    kupr <- kstart
  } else if (!is.null(kbounds)) {
    # validate the bounds, and any starting values, against the data
    if (
      kbounds[1] <= hsort[3] ||
        kbounds[2] >= hsort[n - 2] ||
        any(kstart < kbounds[1]) ||
        any(kstart > kbounds[2])
    ) {
      stop("invalid `kstart` or `kbounds`: need hsort[3] < kbounds[1] <= kstart <= kbounds[2] < hsort[n - 2]")
    }
    klwr <- kbounds[1]
    kupr <- kbounds[2]
  } else {
    klwr <- hsort[3] + 0.001
    kupr <- hsort[n - 2] - 0.001
    if (any(kstart < klwr) || any(kstart > kupr)) {
      stop("`kstart` must leave at least 3 gaugings in each segment")
    }
  }
  if (is.null(kstart)) {
    # 10 interior points of the search range
    kstart <- seq(klwr, kupr, length.out = 12L)[2:11]
  }

  # -- 4. Model formula ----
  # Model formula (see the header for the full derivation). The upper branch
  # of the successive form has no free a2: the leading coefficient
  # a1*(k - c1)^b1 / (k - c2)^b2 is exactly what makes the two branches equal
  # at stage = k, enforcing continuity through the `a` parameter.
  if (controls == "successive") {
    modform <- discharge ~ ifelse(
      stage < k,
      a1 * (stage - c1)^b1,
      (a1 * (k - c1)^b1 / (k - c2)^b2) * (stage - c2)^b2
    )
  } else if (controls == "additive") {
    # Upper branch adds an extra power law to the low-flow discharge at k.
    modform <- discharge ~ ifelse(
      stage < k,
      a1 * (stage - c1)^b1,
      a1 * (k - c1)^b1 + a2 * (stage - k)^b2
    )
  }

  if (wts_code == "none") {
    wts <- rep(1, length(discharge))
  }
  if (wts_code != "prop") {
    checkmate::assert_numeric(wts, len = length(discharge), .var.name = "wts")
  }

  # Steps 3, 5 and 6 depend on the starting breakpoint, so they are wrapped
  # up to be repeated for each one.
  fit_from <- function(kstart) {
    # -- 3. Starting values, from a log-log fit to each segment ----
    # Starting values for the nls fit. Split the data at kstart and fit each
    # segment separately as a straight line on the log-log scale, since
    # log(discharge) = log(a) + b*log(stage - c) is linear in log(a) and b once c is fixed.
    qh1 <- subset(qh, stage < kstart) # low-flow segment
    qh2 <- subset(qh, stage >= kstart) # high-flow segment

    # Lower segment: pick c1 just below the smallest stage so that (stage - c1) > 0,
    # then read a1, b1 off the log-log linear fit.
    c1start <- min(qh1$stage) - 0.1 * (max(qh1$stage) - min(qh1$stage))
    lm_mod <- stats::lm(log(qh1$discharge) ~ log(qh1$stage - c1start))
    pars_1 <- as.numeric(stats::coefficients(lm_mod))
    a1start <- exp(pars_1[1])
    b1start <- pars_1[2]

    # Upper segment: starting values depend on how the segments are joined.
    if (controls == "successive") {
      # Independent power law on the upper data; offset c2 placed between c1 and k.
      c2start <- 0.5 * (kstart + c1start)
      lm_mod <- stats::lm(log(qh2$discharge) ~ log(qh2$stage - c2start))
      pars_2 <- as.numeric(stats::coefficients(lm_mod))
      a2start <- exp(pars_2[1])
      b2start <- pars_2[2]
    } else if (controls == "additive") {
      # Remove the low-flow discharge carried up to the breakpoint, then fit the
      # remaining "excess" discharge against depth above k, (stage - k). Keep only
      # positive residuals so the log is defined.
      if (kstart < c1start) {
        stop("`kstart` is below the starting value of `c1`")
      }
      qh2$excess <- qh2$discharge - a1start * (kstart - c1start)^b1start
      qh2 <- qh2[which(qh2$excess > 0), , drop = FALSE]
      lm_mod <- stats::lm(log(excess) ~ log(stage - kstart), data = qh2)
      pars_2 <- as.numeric(stats::coefficients(lm_mod))
      a2start <- exp(pars_2[1])
      b2start <- pars_2[2]
    }

    # -- 5. Assemble start values and bounds for the port algorithm ----
    # starting values and bounds for nls arguments
    if (controls == "successive") {
      start_list <- list(
        a1 = a1start,
        b1 = b1start,
        c1 = c1start,
        b2 = b2start,
        c2 = c2start,
        k = kstart
      )
      lwr_list <- list(
        a1 = 0,
        b1 = 0,
        c1 = -Inf,
        b2 = 0,
        c2 = -Inf,
        k = klwr
      )
      upr_list <- list(
        a1 = Inf,
        b1 = 4,
        c1 = min(stage) - 0.001,
        b2 = 4,
        c2 = max(stage),
        k = kupr
      )
    } else if (controls == "additive") {
      start_list <- list(
        a1 = a1start,
        b1 = b1start,
        c1 = c1start,
        a2 = a2start,
        b2 = b2start,
        k = kstart
      )
      lwr_list <- list(
        a1 = 0,
        b1 = 0,
        c1 = -Inf,
        a2 = 0,
        b2 = 0,
        k = klwr
      )
      upr_list <- list(
        a1 = Inf,
        b1 = 4,
        c1 = min(stage) - 0.001,
        a2 = Inf,
        b2 = 4,
        k = kupr
      )
    }
    # -- 6. Fit: one pass for fixed weights, IRLS for "prop" ----
    # Fit the model. Two paths: a single nls fit when the weights are known up
    # front ("none" or user-supplied "spec"), or an iteratively reweighted loop
    # when the weights depend on the fitted values ("prop").
    #
    # Path 1: fixed weights (OLS if wts_code == "none", else the supplied wts).
    if (wts_code != "prop") {
      mod_nls <- stats::nls(
        formula = modform,
        data = data.frame(discharge, stage),
        weights = wts,
        start = start_list,
        lower = unlist(lwr_list),
        upper = unlist(upr_list),
        control = control,
        algorithm = "port"
      )
      # Path 2: proportional weights (constant coefficient of variation). The
      # weights 1/fitted^2 depend on the (unknown) fitted discharge, so we iterate:
      # fit -> recompute weights from the new fitted values -> refit, stopping when
      # the fitted discharges change by less than tol (or after maxiter). The
      # initial weights use the log-log starting-value curve.
    } else if (wts_code == "prop") {
      if (controls == "successive") {
        yp <- ifelse(
          stage < kstart,
          a1start * (stage - c1start)^b1start,
          a2start * (stage - c2start)^b2start
        )
      } else {
        # controls = "additive"
        yp <- ifelse(
          stage < kstart,
          a1start * (stage - c1start)^b1start,
          a1start * (kstart - c1start)^b1start + a2start * (stage - kstart)^b2start
        )
      }
      fit_fun <- function(wts, start) {
        stats::nls(
          formula = modform,
          data = data.frame(discharge, stage, wts),
          weights = wts,
          start = start,
          lower = unlist(lwr_list),
          upper = unlist(upr_list),
          control = control,
          algorithm = "port"
        )
      }
      res <- reweight_in_rounds(
        fit_fun,
        yp = yp,
        start = start_list,
        tol = weighting$tol,
        maxiter = weighting$maxiter
      )
      return(res)
    }
    list(model = mod_nls, weights = wts, irls = NULL)
  }

  # -- 7. Try each starting breakpoint, and keep the best fit ----
  # Warnings from individual starts (unconverged reweighting, say) are
  # dropped here; the one that matters, for the fit kept, is reissued below.
  fits <- lapply(kstart, function(ks) {
    tryCatch(suppressWarnings(fit_from(ks)), error = function(e) e)
  })
  ok <- !vapply(fits, inherits, logical(1), what = "error")
  if (!any(ok)) {
    stop(
      "The two-segment fit failed from every starting breakpoint tried. ",
      "The first error was: ", conditionMessage(fits[[1]]),
      call. = FALSE
    )
  }
  loglik <- rep(NA_real_, length(fits))
  k_hat <- rep(NA_real_, length(fits))
  for (i in which(ok)) {
    mu <- as.numeric(stats::fitted(fits[[i]]$model))
    w_model <- if (wts_code == "prop") 1 / mu^2 else wts
    loglik[i] <- profile_loglik(discharge, mu, w_model)
    k_hat[i] <- stats::coef(fits[[i]]$model)[["k"]]
  }
  kstart_search <- tibble::tibble(kstart = kstart, k = k_hat, loglik = loglik)
  # prefer fits whose reweighting converged, then the highest likelihood
  converged <- vapply(
    fits,
    function(f) !inherits(f, "error") && !isFALSE(f$irls$converged),
    logical(1)
  )
  candidates <- if (any(converged)) which(converged) else which(ok)
  best <- candidates[which.max(loglik[candidates])]
  mod_nls <- fits[[best]]$model
  wts <- fits[[best]]$weights
  irls <- fits[[best]]$irls
  if (isFALSE(irls$converged)) {
    warning(
      "Proportional weights did not converge in ", weighting$maxiter, " rounds ",
      "from any starting breakpoint; the fit may not be reliable. ",
      "Consider increasing `maxiter` in `wts_prop()`.",
      call. = FALSE
    )
  }
  qh <- tibble::as_tibble(qh)
  mod_sum <- summary(mod_nls)
  th <- stats::coef(mod_nls)
  pars <- if (controls == "successive") {
    list(
      a = th[["a1"]],
      b = unname(th[c("b1", "b2")]),
      c = unname(th[c("c1", "c2")]),
      k = th[["k"]]
    )
  } else {
    list(
      a = unname(th[c("a1", "a2")]),
      b = unname(th[c("b1", "b2")]),
      c = th[["c1"]],
      k = th[["k"]]
    )
  }
  outlist <- list(
    gaugings = qh,
    pars = pars,
    # Everything needed to refit these data from scratch. boot_limits_2seg()
    # resamples and refits, so it has to reproduce the original call exactly;
    # without this it would silently fall back on the argument defaults.
    settings = list(
      controls = controls,
      kstart = kstart_input,
      kfixed = kfixed,
      kbounds = kbounds,
      wts = weighting,
      control = control
    ),
    weights = wts,
    irls = irls,
    wts = weighting,
    kstart_search = kstart_search,
    rse = mod_sum$sigma,
    model = mod_nls
  )
  structure(outlist, class = c("rc_2seg_power", "rating_curve"))
}
