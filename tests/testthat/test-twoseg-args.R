# Two-segment fitting arguments, and the interval methods under each
# weighting. Sauze (from RBaM) is used where the fit needs a clear change of
# control; the Thompson is enough for argument checking.

sauze <- function() {
  skip_if_not_installed("RBaM")
  RBaM::SauzeGaugings
}

# A piecewise power law with a break at stage 1.5, and deterministic scatter.
two_control <- function() {
  stage <- seq(0.2, 3, length.out = 40)
  k <- 1.5
  mu <- ifelse(
    stage < k,
    10 * stage^1.5,
    10 * k^1.5 / (k - 1)^2.5 * (stage - 1)^2.5
  )
  data.frame(stage = stage, discharge = mu * (1 + 0.03 * sin(7 * seq_along(stage))))
}

sauze_fit <- function(combine = "replace", variance = "none") {
  d <- sauze()
  if (identical(variance, "spec")) {
    variance <- var_spec(d$uQ^2)
  }
  rc_2seg_powerlaw(Q, H, data = d, combine = combine, variance = variance, kstart = 1)
}

# -- arguments -----------------------------------------------------------------

test_that("too few gaugings is an error", {
  expect_error(rc_2seg_powerlaw(1:6, 1:6), "too few")
})

test_that("a misspelled argument is an error, not ignored", {
  expect_error(
    rc_2seg_powerlaw(discharge, stage, data = thompson, kstrat = 2),
    "must be empty"
  )
  fit <- rc_2seg_powerlaw(Q, H, data = sauze(), kstart = 1)
  expect_error(boot_limits_2seg(fit, new_stage = 2, N = 10), "must be empty")
})

test_that("kfixed holds the breakpoint at kstart", {
  expect_error(rc_2seg_powerlaw(Q, H, data = sauze(), kfixed = TRUE), "kstart")
  fit <- rc_2seg_powerlaw(Q, H, data = sauze(), kfixed = TRUE, kstart = 1.8)
  expect_equal(fit$curve_parameters[["k"]], 1.8)
})

test_that("kbounds limit the breakpoint search", {
  # Synthetic gaugings with a clear break at stage 1.5: the Sauze fits are too
  # sensitive to their starting values to test the bounds on reliably across
  # platforms. This fit converges from any start between about 1.05 and 1.75.
  d <- two_control()
  fit <- rc_2seg_powerlaw(discharge, stage, data = d, kbounds = c(1, 2), kstart = 1.5)
  expect_equal(fit$curve_parameters[["k"]], 1.5, tolerance = 0.02)
  fit2 <- rc_2seg_powerlaw(discharge, stage, data = d, kstart = 1.3, kbounds = c(1, 2))
  expect_equal(fit2$settings$kbounds, c(1, 2))
  expect_equal(fit2$curve_parameters[["k"]], fit$curve_parameters[["k"]], tolerance = 1e-4)
  expect_error(
    rc_2seg_powerlaw(discharge, stage, data = d, kstart = 3, kbounds = c(1, 2)),
    "invalid"
  )
})

test_that("kstart must leave three gaugings in each segment", {
  expect_error(
    rc_2seg_powerlaw(Q, H, data = sauze(), kstart = min(sauze()$H)),
    "at least 3"
  )
})

test_that("two-segment fits keep the specified weights for the bootstrap", {
  d <- sauze()
  fit <- rc_2seg_powerlaw(Q, H, data = d, variance = var_spec(d$uQ^2), kstart = 1)
  expect_equal(fit$settings$variance$values, d$uQ^2)
})

test_that("specified variances stay aligned when gaugings are dropped", {
  d <- sauze()
  w <- d$uQ^2
  d$Q[5] <- NA
  fit <- rc_2seg_powerlaw(Q, H, data = d, variance = var_spec(w), kstart = 1)
  expect_equal(fit$settings$variance$values, w[-5])
  # and the bootstrap, which resamples them, still runs
  b <- suppressWarnings(
    boot_limits_2seg(fit, new_stage = 2, conflev = 0.9, B = 5, seed = 1)
  )
  expect_equal(nrow(b), 1L)
})

# -- delta ---------------------------------------------------------------------

test_that("delta prediction limits under each weighting", {
  for (wc in c("none", "prop")) {
    fit <- sauze_fit("replace", wc)
    p <- predict(fit, new_stage = c(1, 3), conflev = 0.95, predlev = 0.95)
    expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr), info = wc)
  }
})

# -- boot ----------------------------------------------------------------------

test_that("boot prediction limits under each weighting", {
  for (wc in c("none", "prop", "spec")) {
    fit <- sauze_fit("replace", wc)
    p <- suppressWarnings(
      boot_limits_2seg(fit, new_stage = c(1, 3), conflev = 0.9, predlev = 0.9,
                       B = 10, seed = 1)
    )
    expect_named(p, c("stage", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr"))
    if (wc == "spec") {
      expect_true(all(is.na(p$pi_lwr)))
    } else {
      expect_true(all(p$pi_lwr < p$fit & p$fit < p$pi_upr), info = wc)
    }
    expect_lte(attr(p, "B_success"), 10L)
  }
})

test_that("boot refuses spec weights it cannot resample", {
  fit <- sauze_fit("replace", "spec")
  fit$settings$variance$values <- NULL
  expect_error(boot_limits_2seg(fit, new_stage = 2, B = 5), "variances given")
})

test_that("boot warns when too few resamples converge", {
  fit <- sauze_fit("replace", "none")
  # corrupt the arguments so that every refit fails
  fit$settings$kstart <- -100
  expect_warning(
    p <- boot_limits_2seg(fit, new_stage = 2, conflev = 0.9, B = 3,
                          max_tries_factor = 1),
    "only 0 of 3"
  )
  expect_equal(attr(p, "B_success"), 0L)
})

# -- starting breakpoints --------------------------------------------------------

test_that("by default, starting breakpoints span the search range", {
  d <- two_control()
  fit <- rc_2seg_powerlaw(discharge, stage, data = d)
  s <- fit$kstart_search
  expect_named(s, c("kstart", "k", "loss"))
  expect_equal(nrow(s), 10L)
  hs <- sort(d$stage)
  expect_true(all(s$kstart > hs[3] & s$kstart < hs[length(hs) - 2]))
  expect_null(fit$settings$kstart)
  # the fit kept has the smallest loss of those that succeeded
  expect_equal(fit$curve_parameters[["k"]], s$k[which.min(s$loss)])
  expect_equal(fit$curve_parameters[["k"]], 1.5, tolerance = 0.02)
})

test_that("kbounds set the range the default starts span", {
  fit <- rc_2seg_powerlaw(discharge, stage, data = two_control(), kbounds = c(1, 2))
  expect_true(all(fit$kstart_search$kstart > 1 & fit$kstart_search$kstart < 2))
  expect_equal(fit$curve_parameters[["k"]], 1.5, tolerance = 0.02)
})

test_that("a vector of starting breakpoints is tried in turn", {
  fit <- rc_2seg_powerlaw(discharge, stage, data = two_control(), kstart = c(1.2, 1.4, 1.6))
  expect_equal(fit$kstart_search$kstart, c(1.2, 1.4, 1.6))
  expect_equal(fit$settings$kstart, c(1.2, 1.4, 1.6))
  expect_error(
    rc_2seg_powerlaw(discharge, stage, data = two_control(), kfixed = TRUE, kstart = c(1, 2)),
    "single value"
  )
})

test_that("starts that fail are recorded and passed over", {
  skip_if_not_installed("RBaM")
  d <- RBaM::SauzeGaugings
  # from starting breakpoints of 1.95 and above this fit used to fail
  fit <- rc_2seg_powerlaw(Q, H, data = d, variance = var_spec(d$uQ^2),
                     kbounds = c(1.5, 2.5))
  expect_true(anyNA(fit$kstart_search$k))
  expect_equal(fit$curve_parameters[["k"]], 1.62, tolerance = 0.01)
})

test_that("a fit that fails from every start is an error", {
  expect_error(
    rc_2seg_powerlaw(-thompson$discharge, thompson$stage),
    "failed from every starting breakpoint"
  )
})

test_that("fits are compared on the loss they minimise", {
  discharge <- c(1, 2, 3, 4)
  # fixed weights: the weighted residual sum of squares
  expect_equal(
    fit_loss(discharge, discharge + 0.1, w = rep(2, 4), variance_type = "none"),
    2 * 4 * 0.01
  )
  # proportional: the negative Gamma quasi-likelihood, smallest at mu = y
  at_data <- fit_loss(discharge, discharge, w = NULL, variance_type = "prop")
  expect_equal(at_data, sum(1 + log(discharge)))
  expect_lt(at_data, fit_loss(discharge, discharge * 1.1, NULL, "prop"))
  expect_lt(at_data, fit_loss(discharge, discharge * 0.9, NULL, "prop"))
})


test_that("the add form's upper starting values never fail", {
  # a gauging exactly at the starting breakpoint is not logged
  qh2 <- tibble::tibble(stage = c(2, 2.5, 3), discharge = c(50, 90, 150))
  start <- add_upper_start(qh2, kstart = 2, a1 = 10, b1 = 1.5, c1 = 0)
  expect_true(all(is.finite(unlist(start))))
  # fewer than two gaugings with a positive excess: start from the lower
  # power law instead of fitting a line
  low <- tibble::tibble(stage = c(2.2, 2.5, 3), discharge = c(10, 20, 80))
  expect_equal(
    add_upper_start(low, kstart = 2, a1 = 10, b1 = 1.5, c1 = 0),
    list(a2 = 10, b2 = 1.5)
  )
})

test_that("a starting breakpoint at a gauged stage works under add", {
  d <- sauze()
  at_gauging <- sort(d$H)[length(d$H) - 3]
  fit <- rc_2seg_powerlaw(Q, H, data = d, combine = "add", kstart = at_gauging)
  expect_true(is.finite(fit$curve_parameters$k))
})
