# Weighting schemes, and the settings passed to nls().

test_that("the shorthand strings match the constructors", {
  a <- rc_power(discharge, stage, data = thompson, wts = "prop")
  b <- rc_power(discharge, stage, data = thompson, wts = wts_prop())
  expect_equal(coef(a), coef(b))
  expect_equal(
    coef(rc_power(discharge, stage, data = thompson, wts = "none")),
    coef(rc_power(discharge, stage, data = thompson))
  )
})

test_that("invalid weighting schemes are errors", {
  expect_error(
    rc_power(discharge, stage, data = thompson, wts = "spec"),
    "wts_spec"
  )
  expect_error(rc_power(discharge, stage, data = thompson, wts = "bogus"))
  expect_error(rc_power(discharge, stage, data = thompson, wts = 1:3), "must be")
  expect_error(wts_prop(tol = -1))
  expect_error(wts_prop(maxiter = 0))
})

test_that("the fit records its weighting scheme", {
  fit <- rc_power(discharge, stage, data = thompson, wts = wts_prop(tol = 1e-4))
  expect_s3_class(fit$settings$wts, "rc_wts_prop")
  expect_equal(fit$settings$wts$tol, 1e-4)
  expect_equal(fit$settings$wts$maxiter, 100)
})

test_that("weighting schemes print what they are", {
  expect_output(print(wts_none()), "none")
  expect_output(print(wts_prop(maxiter = 20)), "maxiter = 20")
  expect_output(print(wts_spec(1 / uncertainty_sd^2)), "1/uncertainty_sd\\^2")
  d <- thompson[!is.na(thompson$uncertainty_pct), ]
  d$uncertainty_sd <- d$uncertainty_pct / 100 * d$discharge / 2
  fit <- rc_power(discharge, stage, data = d,
                wts = wts_spec(1 / uncertainty_sd^2))
  expect_output(print(fit$settings$wts), "for 19 gaugings")
})

test_that("loess reweights proportionally in rounds, like the other fits", {
  fit <- rc_loess(discharge, stage, data = thompson, wts = "prop")
  expect_true(fit$irls$converged)
  expect_gt(fit$irls$iterations, 1L)
  expect_null(rc_loess(discharge, stage, data = thompson)$irls)
  expect_warning(
    rc_loess(discharge, stage, data = thompson, wts = wts_prop(maxiter = 1)),
    "did not converge"
  )
})

test_that("control is passed to nls()", {
  tight <- stats::nls.control(maxiter = 1)
  expect_error(
    rc_power(discharge, stage, data = thompson, control = tight),
    "maximum of 1"
  )
  expect_error(
    rc_power_log(discharge, stage, data = thompson, control = tight),
    "maximum of 1"
  )
  expect_error(
    rc_2seg_power(discharge, stage, data = thompson, wts = "prop",
                control = list(maxiter = 1)),
    "failed from every starting breakpoint"
  )
  expect_error(
    rc_power(discharge, stage, data = thompson, control = list(1)),
    "names"
  )
})

test_that("the old arguments are gone", {
  for (old in list(
    list(wts_code = "prop"), list(wts_tol = 1e-3), list(wts_maxiter = 5),
    list(nls_tol = 1e-3), list(nls_maxiter = 5), list(config = "additive")
  )) {
    expect_error(
      do.call(rc_2seg_power, c(list(thompson$discharge, thompson$stage), old)),
      class = "rlib_error_dots_nonempty"
    )
  }
})

# -- exponent of proportional scatter --------------------------------------------

test_that("an exponent of 1 is the default proportional scatter", {
  expect_equal(
    coef(rc_power(discharge, stage, data = thompson, wts = wts_prop(exponent = 1))),
    coef(rc_power(discharge, stage, data = thompson, wts = "prop"))
  )
})

test_that("a fixed exponent works in every fitting function", {
  w <- wts_prop(exponent = 1.5)
  fits <- list(
    rc_power(discharge, stage, data = thompson, wts = w),
    rc_poly(discharge, stage, data = thompson, wts = w),
    rc_loess(discharge, stage, data = thompson, wts = w),
    rc_2seg_power(discharge, stage, data = thompson, wts = w, kstart = 2)
  )
  for (fit in fits) {
    expect_equal(fit$exponent, 1.5)
    expect_true(fit$irls$converged)
    mu <- as.numeric(stats::predict(fit$model))
    expect_equal(fit$weights, 1 / mu^3, tolerance = 1e-3)
  }
  # with a larger exponent the prediction band widens faster with the flow
  p1 <- predict(rc_power(discharge, stage, data = thompson, wts = "prop"),
                stage = c(1, 6), predlev = 0.95)
  p15 <- predict(fits[[1]], stage = c(1, 6), predlev = 0.95)
  ratio <- function(p) diff(p$pi_upr - p$pi_lwr) / (p$pi_upr - p$pi_lwr)[1]
  expect_gt(ratio(p15), ratio(p1))
})

test_that("the exponent can be estimated in rc_power() only", {
  fit <- rc_power(discharge, stage, data = thompson,
                  wts = wts_prop(exponent = NULL))
  # the same estimate as nlme::gnls() from the same start
  d <- as.data.frame(thompson[, c("discharge", "stage")])
  ref <- nlme::gnls(
    discharge ~ a * (stage - c)^b,
    data = d,
    start = as.list(coef(rc_power(discharge, stage, data = thompson,
                                  wts = "prop"))),
    weights = nlme::varPower(),
    control = nlme::gnlsControl(maxIter = 1e5, minScale = 1e-5)
  )
  expect_equal(coef(fit), coef(ref))
  p <- predict(fit, stage = c(1, 3), conflev = 0.95, predlev = 0.95)
  expect_named(p, c("stage", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr"))
  expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr))
  est <- wts_prop(exponent = NULL)
  for (fitter in list(rc_poly, rc_loess)) {
    expect_error(
      fitter(discharge, stage, data = thompson, wts = est),
      "available in `rc_power\\(\\)` only"
    )
  }
  expect_error(
    rc_2seg_power(discharge, stage, data = thompson, wts = est),
    "available in `rc_power\\(\\)` only"
  )
})

test_that("the exponent is validated and printed", {
  expect_error(wts_prop(exponent = "a"))
  expect_output(print(wts_prop(exponent = NULL)), "exponent estimated")
  expect_output(print(wts_prop(exponent = 2)), "to the power 2")
  expect_output(print(wts_prop()), "proportional to the flow \\(tol")
})
