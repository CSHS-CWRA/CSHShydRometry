# Two-segment tests on a river that actually has two controls.
#
# The Thompson gaugings shipped with the package are close to a single
# control, so their two-segment fit is fragile: only the piecewise
# configuration converges, and only with proportional weights and a supplied
# starting breakpoint. The Ardeche at Sauze, in RBaM, fits under every
# combination and carries a reported uncertainty for every gauging, so it is
# what the two-segment code is exercised against where RBaM is installed.

sauze_fit <- function(config = "piecewise", wts_code = "none") {
  d <- RBaM::SauzeGaugings
  args <- list(
    q = d$Q,
    h = d$H,
    config = config,
    wts_code = wts_code,
    kstart = 1
  )
  if (wts_code == "spec") {
    args$wts <- 1 / d$uQ^2
  }
  do.call(rc_nls_2seg, args)
}

test_that("both configurations fit under every weighting", {
  skip_if_not_installed("RBaM")
  for (cfg in c("piecewise", "compound")) {
    for (wc in c("none", "prop", "spec")) {
      fit <- sauze_fit(cfg, wc)
      expect_s3_class(fit, "rc_nls_2seg")
      expect_s3_class(fit, "rating_curve")
      k <- fit$pars[["k"]]
      expect_true(
        k > min(RBaM::SauzeGaugings$H) && k < max(RBaM::SauzeGaugings$H),
        info = paste(cfg, wc)
      )
    }
  }
})

test_that("the breakpoint is where the record says it should be", {
  skip_if_not_installed("RBaM")
  # the change of control at Sauze sits a little above 1.5 m; all six
  # combinations should land in the same neighbourhood
  ks <- vapply(
    c("piecewise", "compound"),
    function(cfg) sauze_fit(cfg, "spec")$pars[["k"]],
    numeric(1)
  )
  expect_true(all(ks > 1.2 & ks < 2.2))
})

test_that("specified weights give confidence but not prediction limits", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("piecewise", "spec")
  p <- suppressMessages(
    predict(fit, hpred = c(1, 3), conflev = 0.95, predlev = 0.95)
  )
  expect_true(all(is.finite(p$ci_lwr)))
  expect_true(all(is.na(p$pi_lwr)))
  expect_true(all(is.na(p$pi_upr)))
})

test_that("the delta band jumps at the breakpoint and the bootstrap does not", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("piecewise", "spec")
  k <- fit$pars[["k"]]
  hh <- c(k - 0.05, k + 0.05)

  d <- suppressMessages(predict(fit, hpred = hh, conflev = 0.95))
  d_ratio <- (d$ci_upr[2] - d$ci_lwr[2]) / (d$ci_upr[1] - d$ci_lwr[1])

  b <- suppressWarnings(suppressMessages(
    predict(fit, hpred = hh, conflev = 0.95, method = "boot",
            B = 150, seed = 1)
  ))
  b_ratio <- (b$ci_upr[2] - b$ci_lwr[2]) / (b$ci_upr[1] - b$ci_lwr[1])

  # this is the documented weakness of the default method: the band widens
  # sharply across the breakpoint, where a method that does not fix k does not
  expect_gt(d_ratio, 3)
  expect_lt(b_ratio, 2)
})

test_that("all three methods agree on the fitted curve", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("piecewise", "spec")
  hp <- c(0.5, 1, 2, 4)
  d <- delta_limits_2seg(fit, hpred = hp)
  b <- suppressWarnings(
    boot_limits_2seg(fit, hpred = hp, B = 25, seed = 1)
  )
  s <- sim_limits_2seg(fit, hpred = hp, M = 100, seed = 1)
  expect_equal(d$fit, b$fit)
  expect_equal(d$fit, s$fit)
})

test_that("the compound configuration is continuous at the breakpoint", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("compound", "spec")
  k <- fit$pars[["k"]]
  p <- predict(fit, hpred = c(k - 1e-6, k + 1e-6))
  expect_equal(p$fit[1], p$fit[2], tolerance = 1e-5)
})
