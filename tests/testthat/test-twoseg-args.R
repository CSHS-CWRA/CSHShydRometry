# Two-segment fitting arguments, and the interval methods under each
# weighting. Sauze (from RBaM) is used where the fit needs a clear change of
# control; the Thompson is enough for argument checking.

sauze <- function() {
  skip_if_not_installed("RBaM")
  RBaM::SauzeGaugings
}

sauze_fit <- function(config = "piecewise", wts_code = "none") {
  d <- sauze()
  if (wts_code == "spec") {
    rc_nls_2seg(Q, H, data = d, config = config, wts_code = "spec",
                wts = 1 / d$uQ^2, kstart = 1)
  } else {
    rc_nls_2seg(Q, H, data = d, config = config, wts_code = wts_code,
                kstart = 1)
  }
}

# -- arguments -----------------------------------------------------------------

test_that("too few gaugings is an error", {
  expect_error(rc_nls_2seg(1:6, 1:6), "too few")
})

test_that("a misspelled argument is an error, not ignored", {
  expect_error(
    rc_nls_2seg(q, h, data = thompson, kstrat = 2),
    "must be empty"
  )
  fit <- rc_nls_2seg(Q, H, data = sauze(), kstart = 1)
  expect_error(boot_limits_2seg(fit, hpred = 2, N = 10), "must be empty")
})

test_that("the c continuity constraint is not implemented", {
  expect_error(
    rc_nls_2seg(q, h, data = thompson, contcons = "c", kstart = 2),
    "not implemented"
  )
})

test_that("kfixed holds the breakpoint at kstart", {
  expect_error(rc_nls_2seg(Q, H, data = sauze(), kfixed = TRUE), "kstart")
  fit <- rc_nls_2seg(Q, H, data = sauze(), kfixed = TRUE, kstart = 1.8)
  expect_equal(fit$pars[["k"]], 1.8)
})

test_that("kbounds limit the breakpoint search", {
  d <- sauze()
  fit <- rc_nls_2seg(Q, H, data = d, kbounds = c(1.5, 2.5))
  expect_true(fit$pars[["k"]] >= 1.5 && fit$pars[["k"]] <= 2.5)
  fit2 <- rc_nls_2seg(Q, H, data = d, kstart = 1.8, kbounds = c(1.5, 2.5))
  expect_equal(fit2$settings$kbounds, c(1.5, 2.5))
  expect_error(
    rc_nls_2seg(Q, H, data = d, kstart = 3, kbounds = c(1.5, 2.5)),
    "invalid"
  )
})

test_that("kstart must leave three gaugings in each segment", {
  expect_error(
    rc_nls_2seg(Q, H, data = sauze(), kstart = min(sauze()$H)),
    "at least 3"
  )
})

test_that("specified weights stay aligned when gaugings are dropped", {
  d <- sauze()
  w <- 1 / d$uQ^2
  d$Q[5] <- NA
  fit <- rc_nls_2seg(Q, H, data = d, wts_code = "spec", wts = w, kstart = 1)
  expect_equal(fit$settings$wts, w[-5])
  # and the bootstrap, which resamples them, still runs
  b <- suppressWarnings(
    boot_limits_2seg(fit, hpred = 2, conflev = 0.9, B = 5, seed = 1)
  )
  expect_equal(nrow(b), 1L)
})

# -- delta ---------------------------------------------------------------------

test_that("delta prediction limits under each weighting", {
  for (wc in c("none", "prop")) {
    fit <- sauze_fit("piecewise", wc)
    p <- predict(fit, hpred = c(1, 3), conflev = 0.95, predlev = 0.95)
    expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr), info = wc)
  }
})

# -- boot ----------------------------------------------------------------------

test_that("boot prediction limits under each weighting", {
  for (wc in c("none", "prop", "spec")) {
    fit <- sauze_fit("piecewise", wc)
    p <- suppressWarnings(
      boot_limits_2seg(fit, hpred = c(1, 3), conflev = 0.9, predlev = 0.9,
                       B = 10, seed = 1)
    )
    expect_named(p, c("h", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr"))
    if (wc == "spec") {
      expect_true(all(is.na(p$pi_lwr)))
    } else {
      expect_true(all(p$pi_lwr < p$fit & p$fit < p$pi_upr), info = wc)
    }
    expect_lte(attr(p, "B_success"), 10L)
  }
})

test_that("boot refuses spec weights it cannot resample", {
  fit <- sauze_fit("piecewise", "spec")
  fit$settings$wts <- NULL
  expect_error(boot_limits_2seg(fit, hpred = 2, B = 5), "spec weights")
})

test_that("boot warns when too few resamples converge", {
  fit <- sauze_fit("piecewise", "none")
  # corrupt the arguments so that every refit fails
  fit$settings$kstart <- -100
  expect_warning(
    p <- boot_limits_2seg(fit, hpred = 2, conflev = 0.9, B = 3,
                          max_tries_factor = 1),
    "only 0 of 3"
  )
  expect_equal(attr(p, "B_success"), 0L)
})
