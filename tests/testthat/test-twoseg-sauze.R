# Two-segment tests on a river that actually has two controls.
#
# The Thompson gaugings shipped with the package are close to a single
# control, so their two-segment fit is fragile: it converges only with
# proportional weights. The Ardeche at Sauze, in RBaM, fits under every
# combination and carries a reported uncertainty for every gauging, so it is
# what the two-segment code is exercised against where RBaM is installed.

sauze_fit <- function(combine = "replace", variance = "none") {
  d <- RBaM::SauzeGaugings
  if (identical(variance, "spec")) {
    variance <- var_spec(d$uQ^2)
  }
  rc_2seg_powerlaw(Q, H, data = d, combine = combine, variance = variance, kstart = 1)
}

test_that("both ways of combining the segments fit under every weighting", {
  skip_if_not_installed("RBaM")
  for (cfg in c("replace", "add")) {
    for (wc in c("none", "prop", "spec")) {
      fit <- sauze_fit(cfg, wc)
      expect_s3_class(fit, "rc_2seg_powerlaw")
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
    c("replace", "add"),
    function(cfg) sauze_fit(cfg, "spec")$pars[["k"]],
    numeric(1)
  )
  expect_true(all(ks > 1.2 & ks < 2.2))
})

test_that("specified weights give confidence but not prediction limits", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("replace", "spec")
  p <- suppressMessages(
    predict(fit, new_stage = c(1, 3), conflev = 0.95, predlev = 0.95)
  )
  expect_true(all(is.finite(p$ci_lwr)))
  expect_true(all(is.na(p$pi_lwr)))
  expect_true(all(is.na(p$pi_upr)))
})

test_that("the delta band jumps at the breakpoint and the bootstrap does not", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("replace", "spec")
  k <- fit$pars[["k"]]
  hh <- c(k - 0.05, k + 0.05)

  d <- suppressMessages(predict(fit, new_stage = hh, conflev = 0.95))
  d_ratio <- (d$ci_upr[2] - d$ci_lwr[2]) / (d$ci_upr[1] - d$ci_lwr[1])

  b <- suppressWarnings(suppressMessages(
    predict(fit, new_stage = hh, conflev = 0.95, method = "boot",
            B = 150, seed = 1)
  ))
  b_ratio <- (b$ci_upr[2] - b$ci_lwr[2]) / (b$ci_upr[1] - b$ci_lwr[1])

  # this is the documented weakness of the default method: the band widens
  # sharply across the breakpoint, where a method that does not fix k does not
  expect_gt(d_ratio, 3)
  expect_lt(b_ratio, 2)
})

test_that("both methods agree on the fitted curve", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("replace", "spec")
  hp <- c(0.5, 1, 2, 4)
  d <- delta_limits_2seg(fit, new_stage = hp)
  b <- suppressWarnings(
    boot_limits_2seg(fit, new_stage = hp, B = 25, seed = 1)
  )
  expect_equal(d$fit, b$fit)
})

test_that("adding the segments is continuous at the breakpoint", {
  skip_if_not_installed("RBaM")
  fit <- sauze_fit("add", "spec")
  k <- fit$pars[["k"]]
  p <- predict(fit, new_stage = c(k - 1e-6, k + 1e-6))
  expect_equal(p$fit[1], p$fit[2], tolerance = 1e-5)
})
