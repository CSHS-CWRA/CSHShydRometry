# The structure every fit shares: pars, settings, coef(), and the record of
# the proportional-weight reweighting.

all_fits <- function() {
  list(
    rc_log_ols = rc_log_ols(discharge, stage, data = thompson),
    rc_log_nls = rc_log_nls(discharge, stage, data = thompson),
    rc_nls = rc_nls(discharge, stage, data = thompson),
    rc_nls_prop = rc_nls(discharge, stage, data = thompson, wts = "prop"),
    rc_gnls = rc_gnls(discharge, stage, data = thompson),
    rc_poly = rc_poly(discharge, stage, data = thompson),
    rc_poly_prop = rc_poly(discharge, stage, data = thompson, wts = "prop"),
    rc_loess = rc_loess(discharge, stage, data = thompson),
    rc_2seg_nls = rc_2seg_nls(discharge, stage, data = thompson, wts = "prop",
                              kstart = 2)
  )
}

test_that("every fit has the shared elements", {
  for (nm in names(fits <- all_fits())) {
    fit <- fits[[nm]]
    expect_s3_class(fit, "rating_curve")
    expect_named(fit$gaugings, c("discharge", "stage"), info = nm)
    expect_type(fit$pars, "list")
    expect_type(fit$settings, "list")
    expect_true(is.numeric(fit$rse) && length(fit$rse) == 1L, info = nm)
    expect_false(is.null(fit$model), info = nm)
  }
})

test_that("coef() flattens pars into a named numeric vector", {
  for (nm in names(fits <- all_fits())) {
    cf <- coef(fits[[nm]])
    expect_type(cf, "double")
    expect_length(cf, length(unlist(fits[[nm]]$pars)))
    expect_false(is.null(names(cf)), info = nm)
  }
  expect_named(coef(rc_nls(discharge, stage, data = thompson)), c("a", "b", "c"))
  expect_length(coef(rc_loess(discharge, stage, data = thompson)), 0L)
})

test_that("one-segment coef() agrees with the underlying nls model", {
  fit <- rc_nls(discharge, stage, data = thompson)
  expect_equal(coef(fit), stats::coef(fit$model))
})

test_that("two-segment pars hold one value per segment", {
  skip_if_not_installed("RBaM")
  d <- RBaM::SauzeGaugings
  pw <- rc_2seg_nls(Q, H, data = d, kstart = 1)
  expect_named(pw$pars, c("a", "b", "c", "k"))
  expect_equal(lengths(pw$pars), c(a = 1L, b = 2L, c = 2L, k = 1L))
  expect_named(coef(pw), c("a1", "b1", "c1", "b2", "c2", "k"))
  expect_equal(unname(coef(pw)[["b2"]]), pw$pars$b[2])

  cp <- rc_2seg_nls(Q, H, data = d, controls = "additive", kstart = 1)
  expect_equal(lengths(cp$pars), c(a = 2L, b = 2L, c = 1L, k = 1L))
  expect_named(coef(cp), c("a1", "b1", "c1", "a2", "b2", "k"))
})

test_that("settings are enough to refit", {
  fit <- rc_nls(discharge, stage, data = thompson, wts = "prop")
  refit <- do.call(rc_nls, c(list(thompson$discharge, thompson$stage), fit$settings))
  expect_equal(coef(refit), coef(fit))
  fit2 <- rc_poly(discharge, stage, data = thompson, degree = 3)
  refit2 <- do.call(rc_poly, c(list(thompson$discharge, thompson$stage), fit2$settings))
  expect_equal(coef(refit2), coef(fit2))
})

test_that("the two-segment fit no longer takes conflev or predlev", {
  expect_error(
    rc_2seg_nls(discharge, stage, data = thompson, kstart = 2, wts = "prop",
                conflev = 0.9),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("loess records its residual scale and equivalent parameters", {
  fit <- rc_loess(discharge, stage, data = thompson)
  expect_gt(fit$rse, 0)
  expect_gt(fit$enp, 1)
})

test_that("gnls records its variance parameters", {
  fit <- rc_gnls(discharge, stage, data = thompson)
  expect_named(fit$var_pars, "power")
})

# -- reweighting ---------------------------------------------------------------

test_that("converged reweighting is recorded", {
  for (fit in list(
    rc_nls(discharge, stage, data = thompson, wts = "prop"),
    rc_poly(discharge, stage, data = thompson, wts = "prop"),
    rc_2seg_nls(discharge, stage, data = thompson, wts = "prop", kstart = 2)
  )) {
    expect_true(fit$irls$converged)
    expect_gt(fit$irls$iterations, 1L)
    # the weights are those the final model was fitted with
    expect_equal(unname(fit$weights), unname(stats::weights(fit$model)))
  }
  expect_null(rc_nls(discharge, stage, data = thompson)$irls)
})

test_that("reweighting that runs out of rounds warns and says so", {
  expect_warning(
    fit <- rc_nls(discharge, stage, data = thompson, wts = wts_prop(maxiter = 1)),
    "did not converge in 1 rounds"
  )
  expect_false(fit$irls$converged)
  expect_equal(fit$irls$iterations, 1L)
  expect_warning(
    rc_poly(discharge, stage, data = thompson, wts = wts_prop(maxiter = 1)),
    "did not converge"
  )
  expect_warning(
    rc_2seg_nls(discharge, stage, data = thompson, wts = wts_prop(maxiter = 1), kstart = 2),
    "did not converge"
  )
})

test_that("reweighting stops on the change in fitted discharge", {
  fit <- rc_nls(discharge, stage, data = thompson, wts = wts_prop(tol = 1e-3))
  tight <- rc_nls(discharge, stage, data = thompson, wts = wts_prop(tol = 1e-10))
  expect_lte(fit$irls$iterations, tight$irls$iterations)
  expect_equal(stats::fitted(fit$model), stats::fitted(tight$model),
               tolerance = 1e-2)
})
