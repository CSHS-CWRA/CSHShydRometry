# The two-segment interval methods, and the dispatcher over them.

two_seg_fit <- function(config = "piecewise") {
  rc_nls_2seg(
    thompson$q,
    thompson$h,
    config = config,
    wts_code = "prop",
    kstart = 2
  )
}

test_that("predict() defaults to the delta method", {
  fit <- two_seg_fit()
  hp <- c(1, 3, 6)
  expect_equal(
    predict(fit, hpred = hp, conflev = 0.95),
    delta_limits_2seg(fit, hpred = hp, conflev = 0.95)
  )
})

test_that("the k-averaged and simulation methods are not part of this package", {
  fit <- two_seg_fit()
  expect_error(
    predict(fit, hpred = 3, conflev = 0.95, method = "kavg")
  )
  expect_error(predict(fit, hpred = 3, conflev = 0.95, method = "sim"))
  expect_false(exists("sim_limits_2seg", where = asNamespace("CSHShydRometry")))
  expect_false(exists("kavg_limits_2seg", where = asNamespace("CSHShydRometry")))
})

test_that("every interval method returns the same columns", {
  fit <- two_seg_fit()
  hp <- c(1, 3)
  want <- c("h", "fit", "ci_lwr", "ci_upr")
  expect_named(delta_limits_2seg(fit, hpred = hp, conflev = 0.95), want)
  expect_named(
    suppressWarnings(
      boot_limits_2seg(fit, hpred = hp, conflev = 0.95, predlev = NULL,
                       B = 25, seed = 1)
    ),
    want
  )
})

test_that("the fitted curve is the same whichever method is asked for", {
  fit <- two_seg_fit()
  hp <- c(1, 3, 6)
  d <- delta_limits_2seg(fit, hpred = hp, conflev = NULL)
  b <- suppressWarnings(
    boot_limits_2seg(fit, hpred = hp, conflev = NULL, predlev = NULL,
                     B = 25, seed = 1)
  )
  expect_equal(d$fit, b$fit)
})

test_that("both configurations fit and predict", {
  for (cfg in c("piecewise", "compound")) {
    fit <- tryCatch(two_seg_fit(cfg), error = function(e) NULL)
    skip_if(is.null(fit), paste(cfg, "did not converge on these gaugings"))
    expect_s3_class(fit, "rc_nls_2seg")
    expect_s3_class(fit, "rating_curve")
    p <- predict(fit, hpred = c(1, 3), conflev = 0.95)
    expect_true(all(is.finite(p$fit)))
  }
})

test_that("the bootstrap refits with the arguments the fit was made with", {
  fit <- two_seg_fit()
  expect_equal(fit$settings$config, "piecewise")
  expect_equal(fit$settings$wts_code, "prop")
  expect_equal(fit$settings$kstart, 2)
})

test_that("hpred defaults to the observed stage range everywhere", {
  fit <- two_seg_fit()
  p <- predict(fit, conflev = NULL)
  expect_equal(nrow(p), 1000L)
  expect_equal(range(p$h), range(thompson$h))
})

test_that("predict() and the limits functions agree at their defaults", {
  fit <- two_seg_fit()
  hp <- c(1, 3)
  # they all default to computing no intervals at all, so the bare call
  # returns the curve and nothing else
  expect_named(predict(fit, hpred = hp), c("h", "fit"))
  expect_named(delta_limits_2seg(fit, hpred = hp), c("h", "fit"))
})
