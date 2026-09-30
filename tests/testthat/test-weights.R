# Weighting schemes, and the settings passed to nls().

test_that("the shorthand strings match the constructors", {
  a <- rc_nls(discharge, stage, data = thompson, wts = "prop")
  b <- rc_nls(discharge, stage, data = thompson, wts = wts_prop())
  expect_equal(coef(a), coef(b))
  expect_equal(
    coef(rc_nls(discharge, stage, data = thompson, wts = "none")),
    coef(rc_nls(discharge, stage, data = thompson))
  )
})

test_that("invalid weighting schemes are errors", {
  expect_error(
    rc_nls(discharge, stage, data = thompson, wts = "spec"),
    "wts_spec"
  )
  expect_error(rc_nls(discharge, stage, data = thompson, wts = "bogus"))
  expect_error(rc_nls(discharge, stage, data = thompson, wts = 1:3), "must be")
  expect_error(wts_prop(tol = -1))
  expect_error(wts_prop(maxiter = 0))
})

test_that("the fit records its weighting scheme", {
  fit <- rc_nls(discharge, stage, data = thompson, wts = wts_prop(tol = 1e-4))
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
  fit <- rc_nls(discharge, stage, data = d,
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
    rc_nls(discharge, stage, data = thompson, control = tight),
    "maximum of 1"
  )
  expect_error(
    rc_log_nls(discharge, stage, data = thompson, control = tight),
    "maximum of 1"
  )
  expect_error(
    rc_2seg_nls(discharge, stage, data = thompson, wts = "prop",
                control = list(maxiter = 1)),
    "failed from every starting breakpoint"
  )
  expect_error(
    rc_nls(discharge, stage, data = thompson, control = list(1)),
    "names"
  )
})

test_that("the old arguments are gone", {
  for (old in list(
    list(wts_code = "prop"), list(wts_tol = 1e-3), list(wts_maxiter = 5),
    list(nls_tol = 1e-3), list(nls_maxiter = 5), list(config = "additive")
  )) {
    expect_error(
      do.call(rc_2seg_nls, c(list(thompson$discharge, thompson$stage), old)),
      class = "rlib_error_dots_nonempty"
    )
  }
})
