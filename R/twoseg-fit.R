# Fitting the two-segment power-law rating curve.

#' Fit two-segment power-law rating curve using nls on untransformed data
#'
#' @param q Streamflow. A vector, or a column of `data`.
#' @param h Stage. A vector, or a column of `data`.
#' @param data Optional data frame in which to look up `q` and `h`. When
#'   supplied, they may be given as bare column names.
#' @param ... Must be empty. Present so that every argument after it has
#'   to be named, which keeps calls readable and guards against
#'   positional mistakes.
#' @param config Segment configuration, `"piecewise"` or `"compound"`.
#'   Defaults to the first.
#' @param contcons Parameter carrying the continuity constraint, `"a"` or
#'   `"c"`. `"c"` is not implemented.
#' @param kfixed If `TRUE`, hold the breakpoint `k` fixed at `kstart`.
#' @param kstart Starting value for the breakpoint `k`, or `NULL` for the
#'   midpoint of the default search range.
#' @param kbounds Lower and upper bounds for `k`, or `NULL` to keep at least
#'   three gaugings in each segment.
#' @param wts_code Weighting scheme:
#'   `"none"` (ordinary least squares, the default),
#'   `"spec"` (user-supplied weights via `wts`, e.g. `1/uq^2` from reported
#'   discharge uncertainties), or `"prop"` (proportional / constant-CV error,
#'   fit by iteratively reweighting with weights `1/fitted^2`).
#' @param wts Optional vector of weights when `wts_code = "spec"`.
#' @param wts_tol Convergence tolerance under `wts_code = "prop"`: the
#'   reweighting stops once no fitted discharge changes by more than this
#'   fraction from one round to the next.
#' @param wts_maxiter Maximum number of reweighting rounds under
#'   `wts_code = "prop"`. Reaching it gives a warning.
#' @param nls_tol Tolerance for nls convergence.
#' @param nls_maxiter Maximum nls iterations.
#'
#' @details
#' Fitting proceeds in three steps:
#' \enumerate{
#'   \item Choose the breakpoint search range and a starting value for `k`
#'         (from `kstart`/`kbounds`, or defaults that leave >= 3 points per
#'         segment).
#'   \item Derive starting values for the segment parameters by splitting the
#'         data at `kstart` and fitting a linear model to each segment on the
#'         log-log scale (`log q ~ log(h - c)`); `a = exp(intercept)`,
#'         `b = slope`.
#'   \item Fit all parameters jointly with [stats::nls()] using the "port"
#'         algorithm (which supports the box constraints in `lower`/`upper`).
#'         For `wts_code = "prop"` this fit is repeated, updating the weights
#'         from the current fitted values and starting from the previous
#'         estimates, until the fitted discharges stabilise.
#' }
#'
#' @return An object of class `c("rc_nls_2seg", "rating_curve")`; see
#'   [rating_curve] for its contents. `pars` holds the estimated parameters by
#'   type, one value per segment: under `"piecewise"`, `a` has a single value
#'   because the upper segment's coefficient is fixed by continuity, and under
#'   `"compound"`, `c` has a single value because the upper segment is
#'   measured from `k`. `coef()` gives the same estimates by their model
#'   names (`a1`, `b1`, `c1`, ...).
#' @examples
#' # The Thompson is close to a single control, so its two-segment fit needs a
#' # starting breakpoint and proportional weights to converge.
#' fit <- rc_nls_2seg(
#'   q,
#'   h,
#'   data = thompson,
#'   wts_code = "prop",
#'   kstart = 2
#' )
#' fit
#' coef(fit)
#'
#' predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
#'
#' # A river with a clearer change of control is far less fussy. The Ardeche
#' # at Sauze, in the RBaM package, fits under either configuration and
#' # carries a reported uncertainty for every gauging.
#' if (requireNamespace("RBaM", quietly = TRUE)) {
#'   sauze <- RBaM::SauzeGaugings
#'   pw <- rc_nls_2seg(Q, H, data = sauze, kstart = 1)
#'   cp <- rc_nls_2seg(Q, H, data = sauze, config = "compound", kstart = 1)
#'   c(piecewise = pw$pars[["k"]], compound = cp$pars[["k"]])
#'
#'   # weights from the reported gauging uncertainties
#'   rc_nls_2seg(
#'     Q,
#'     H,
#'     data = sauze,
#'     wts_code = "spec",
#'     wts = 1 / sauze$uQ^2,
#'     kstart = 1
#'   )
#' }
#' @export
rc_nls_2seg <- function(
  q,
  h,
  ...,
  data = NULL,
  config = c("piecewise", "compound"),
  contcons = c("a", "c"),
  kfixed = FALSE,
  kstart = NULL,
  kbounds = NULL,
  wts_code = c("none", "spec", "prop"),
  wts = NULL,
  wts_tol = 1e-6,
  wts_maxiter = 100,
  nls_tol = 1e-6,
  nls_maxiter = 1000
) {
  # -- 1. Inputs: tidy evaluation, checks, missing values ----
  # q and h may name columns of `data`, or be vectors
  rlang::check_dots_empty()
  checkmate::assert_data_frame(data, null.ok = TRUE)
  q <- rlang::eval_tidy(rlang::enquo(q), data)
  h <- rlang::eval_tidy(rlang::enquo(h), data)

  # error checks
  config <- rlang::arg_match(config)
  contcons <- rlang::arg_match(contcons)
  wts_code <- rlang::arg_match(wts_code)
  checkmate::assert_numeric(q, min.len = 1L)
  checkmate::assert_numeric(h, len = length(q))
  checkmate::assert_flag(kfixed)
  checkmate::assert_number(kstart, null.ok = TRUE)
  checkmate::assert_numeric(kbounds, len = 2L, null.ok = TRUE)
  checkmate::assert_number(wts_tol, lower = 0)
  checkmate::assert_count(wts_maxiter, positive = TRUE)
  checkmate::assert_number(nls_tol, lower = 0)
  checkmate::assert_count(nls_maxiter, positive = TRUE)
  # remove missing values, keeping user-supplied weights aligned with the
  # gaugings that remain, and check the number of observations
  if (length(wts) == length(q)) {
    wts <- wts[stats::complete.cases(q, h)]
  }
  qh <- rc_complete(q, h)
  # Keep the weights as supplied. `wts` is overwritten below (unit weights for
  # "none", the converged IRLS weights for "prop"), and refitting methods such
  # as boot_limits_2seg need the original to reproduce the fit.
  wts_input <- wts
  irls <- NULL
  q <- qh$q
  h <- qh$h
  hsort <- sort(h)
  n <- length(h)
  if (n < 7) {
    stop("n < 7 - too few data points to fit two-segment curve")
  }
  # -- 2. Breakpoint: search range and starting value ----
  # Determine the breakpoint search range (klwr, kupr) and starting value
  # (kstart). Five cases depending on which of kstart / kbounds the user gave.
  # Defaults keep the breakpoint away from the extreme stages so each segment
  # retains at least 3 observations (hsort[3] .. hsort[n-2]).
  if (kfixed) {
    # hold k fixed at the supplied value (lower == upper bound)
    if (is.null(kstart)) {
      stop("`kstart` must be supplied when `kfixed = TRUE`")
    }
    klwr <- kstart
    kupr <- kstart
  } else if (!is.null(kstart) && !is.null(kbounds)) {
    # both supplied: validate them against the data
    if (
      kbounds[1] <= hsort[3] ||
        kbounds[2] >= hsort[n - 2] ||
        kstart < kbounds[1] ||
        kstart > kbounds[2]
    ) {
      stop("invalid `kstart` or `kbounds`: need hsort[3] < kbounds[1] <= kstart <= kbounds[2] < hsort[n - 2]")
    }
    # use the supplied bounds as the k search range (kstart kept as given)
    klwr <- kbounds[1]
    kupr <- kbounds[2]
  } else if (is.null(kstart) && is.null(kbounds)) {
    # require at least 3 observations for each segment
    klwr <- hsort[3] + 0.001
    kupr <- hsort[n - 2] - 0.001
    kstart <- 0.5 * (klwr + kupr)
  } else if (is.null(kstart) && !is.null(kbounds)) {
    # set kstart to mean of kbounds
    klwr <- kbounds[1]
    kupr <- kbounds[2]
    kstart <- mean(kbounds)
  } else if (!is.null(kstart) && is.null(kbounds)) {
    klwr <- hsort[3] + 0.001
    kupr <- hsort[n - 2] - 0.001
    if (kstart < klwr || kstart > kupr) stop("`kstart` must leave at least 3 gaugings in each segment")
  }
  # -- 3. Starting values, from a log-log fit to each segment ----
  # Starting values for the nls fit. Split the data at kstart and fit each
  # segment separately as a straight line on the log-log scale, since
  # log(q) = log(a) + b*log(h - c) is linear in log(a) and b once c is fixed.
  qh1 <- subset(qh, h < kstart) # low-flow segment
  qh2 <- subset(qh, h >= kstart) # high-flow segment

  # Lower segment: pick c1 just below the smallest stage so that (h - c1) > 0,
  # then read a1, b1 off the log-log linear fit.
  c1start <- min(qh1$h) - 0.1 * (max(qh1$h) - min(qh1$h))
  lm_mod <- stats::lm(log(qh1$q) ~ log(qh1$h - c1start))
  pars_1 <- as.numeric(stats::coefficients(lm_mod))
  a1start <- exp(pars_1[1])
  b1start <- pars_1[2]

  # Upper segment: starting values depend on how the segments are joined.
  if (config == "piecewise") {
    # Independent power law on the upper data; offset c2 placed between c1 and k.
    c2start <- 0.5 * (kstart + c1start)
    lm_mod <- stats::lm(log(qh2$q) ~ log(qh2$h - c2start))
    pars_2 <- as.numeric(stats::coefficients(lm_mod))
    a2start <- exp(pars_2[1])
    b2start <- pars_2[2]
  } else if (config == "compound") {
    # Remove the low-flow discharge carried up to the breakpoint, then fit the
    # remaining "excess" discharge q2 against depth above k, (h - k). Keep only
    # positive residuals so the log is defined.
    if (kstart < c1start) {
      stop("`kstart` is below the starting value of `c1`")
    }
    qh2$q2 <- qh2$q - a1start * (kstart - c1start)^b1start
    qh2 <- qh2[which(qh2$q2 > 0), , drop = FALSE]
    lm_mod <- stats::lm(log(q2) ~ log(h - kstart), data = qh2)
    pars_2 <- as.numeric(stats::coefficients(lm_mod))
    a2start <- exp(pars_2[1])
    b2start <- pars_2[2]
  }

  # -- 4. Model formula ----
  # Model formula (see the header for the full derivation). The upper branch
  # of the piecewise form has no free a2: the leading coefficient
  # a1*(k - c1)^b1 / (k - c2)^b2 is exactly what makes the two branches equal
  # at h = k, enforcing continuity through the `a` parameter (contcons = "a").
  if (config == "piecewise" && contcons == "a") {
    modform <- q ~ ifelse(
      h < k,
      a1 * (h - c1)^b1,
      (a1 * (k - c1)^b1 / (k - c2)^b2) * (h - c2)^b2
    )
  } else if (config == "piecewise" && contcons == "c") {
    stop(
      "`contcons = \"c\"` is not implemented yet"
    )
  } else if (config == "compound") {
    # Upper branch adds an extra power law to the low-flow discharge at k.
    modform <- q ~ ifelse(
      h < k,
      a1 * (h - c1)^b1,
      a1 * (k - c1)^b1 + a2 * (h - k)^b2
    )
  }

  # -- 5. Assemble start values and bounds for the port algorithm ----
  # starting values and bounds for nls arguments
  if (config == "piecewise") {
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
      c1 = min(h) - 0.001,
      b2 = 4,
      c2 = max(h),
      k = kupr
    )
  } else if (config == "compound") {
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
      c1 = min(h) - 0.001,
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
    # create wts vector
    if (wts_code == "none") {
      wts <- rep(1, length(q))
    }
    checkmate::assert_numeric(wts, len = length(q), .var.name = "wts")
    # fit model
    mod_nls <- stats::nls(
      formula = modform,
      data = data.frame(q, h),
      weights = wts,
      start = start_list,
      lower = unlist(lwr_list),
      upper = unlist(upr_list),
      control = list(tol = nls_tol, maxiter = nls_maxiter),
      algorithm = "port"
    )
    # Path 2: proportional weights (constant coefficient of variation). The
    # weights 1/q^2 depend on the (unknown) fitted discharge, so we iterate:
    # fit -> recompute weights from the new fitted values -> refit, stopping when
    # the coefficients change by less than wts_tol (or after wts_maxiter). The
    # initial weights use the log-log starting-value curve.
  } else if (wts_code == "prop") {
    if (config == "piecewise") {
      yp <- ifelse(
        h < kstart,
        a1start * (h - c1start)^b1start,
        a2start * (h - c2start)^b2start
      )
    } else {
      # config = "compound"
      yp <- ifelse(
        h < kstart,
        a1start * (h - c1start)^b1start,
        a1start * (kstart - c1start)^b1start + a2start * (h - kstart)^b2start
      )
    }
    fit_fun <- function(wts, start) {
      stats::nls(
        formula = modform,
        data = data.frame(q, h, wts),
        weights = wts,
        start = start,
        lower = unlist(lwr_list),
        upper = unlist(upr_list),
        control = list(tol = nls_tol, maxiter = nls_maxiter),
        algorithm = "port"
      )
    }
    res <- rc_irls(
      fit_fun,
      yp = yp,
      start = start_list,
      wts_tol = wts_tol,
      wts_maxiter = wts_maxiter
    )
    mod_nls <- res$model
    wts <- res$weights
    irls <- res$irls
  }
  if (requireNamespace("tibble", quietly = TRUE)) {
    qh <- tibble::as_tibble(qh)
  }
  mod_sum <- summary(mod_nls)
  th <- stats::coef(mod_nls)
  pars <- if (config == "piecewise") {
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
    qh_obs = qh,
    pars = pars,
    # Everything needed to refit these data from scratch. boot_limits_2seg()
    # resamples and refits, so it has to reproduce the original call exactly;
    # without this it would silently fall back on the argument defaults.
    settings = list(
      config = config,
      contcons = contcons,
      kfixed = kfixed,
      kstart = kstart,
      kbounds = kbounds,
      wts_code = wts_code,
      wts = wts_input,
      wts_tol = wts_tol,
      wts_maxiter = wts_maxiter,
      nls_tol = nls_tol,
      nls_maxiter = nls_maxiter
    ),
    weights = wts,
    irls = irls,
    rse = mod_sum$sigma,
    model = mod_nls
  )
  structure(outlist, class = c("rc_nls_2seg", "rating_curve"))
}
