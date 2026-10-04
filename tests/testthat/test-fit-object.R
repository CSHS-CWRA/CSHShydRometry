# The structure every fit shares: curve_parameters, settings, coef(), and the record of
# the proportional-weight reweighting.

all_fits <- function() {
  list(
    rc_powerlaw_log_offset = rc_powerlaw_log(discharge, stage, data = thompson, offset = -1.3),
    rc_powerlaw_log = rc_powerlaw_log(discharge, stage, data = thompson),
    rc_powerlaw = rc_powerlaw(discharge, stage, data = thompson),
    rc_power_prop = rc_powerlaw(discharge, stage, data = thompson, variance = "prop"),
    rc_power_var_power = rc_powerlaw(discharge, stage, data = thompson,
                                  variance = var_power()),
    rc_poly = rc_poly(discharge, stage, data = thompson),
    rc_poly_prop = rc_poly(discharge, stage, data = thompson, variance = "prop"),
    rc_loess = rc_loess(discharge, stage, data = thompson),
    rc_2seg_powerlaw = rc_2seg_powerlaw(discharge, stage, data = thompson, variance = "prop",
                              kstart = 2)
  )
}

test_that("every fit has the shared elements", {
  for (nm in names(fits <- all_fits())) {
    fit <- fits[[nm]]
    expect_s3_class(fit, "rating_curve")
    expect_named(fit$gaugings, c("discharge", "stage"), info = nm)
    expect_type(fit$curve_parameters, "list")
    expect_type(fit$settings, "list")
    expect_true(is.numeric(fit$rse) && length(fit$rse) == 1L, info = nm)
    expect_false(is.null(fit$model), info = nm)
  }
})

test_that("coef() flattens curve_parameters into a named numeric vector", {
  for (nm in names(fits <- all_fits())) {
    cf <- coef(fits[[nm]])
    expect_type(cf, "double")
    expect_length(cf, length(unlist(fits[[nm]]$curve_parameters)))
    expect_false(is.null(names(cf)), info = nm)
  }
  expect_named(coef(rc_powerlaw(discharge, stage, data = thompson)), c("a", "b", "c"))
  expect_length(coef(rc_loess(discharge, stage, data = thompson)), 0L)
})

test_that("one-segment coef() agrees with the underlying nls model", {
  fit <- rc_powerlaw(discharge, stage, data = thompson)
  expect_equal(coef(fit), stats::coef(fit$model))
})

test_that("two-segment curve_parameters hold one value per segment", {
  skip_if_not_installed("RBaM")
  d <- RBaM::SauzeGaugings
  pw <- rc_2seg_powerlaw(Q, H, data = d, kstart = 1)
  expect_named(pw$curve_parameters, c("a", "b", "c", "k"))
  expect_equal(lengths(pw$curve_parameters), c(a = 1L, b = 2L, c = 2L, k = 1L))
  expect_named(coef(pw), c("a1", "b1", "c1", "b2", "c2", "k"))
  expect_equal(unname(coef(pw)[["b2"]]), pw$curve_parameters$b[2])

  cp <- rc_2seg_powerlaw(Q, H, data = d, combine = "add", kstart = 1)
  expect_equal(lengths(cp$curve_parameters), c(a = 2L, b = 2L, c = 1L, k = 1L))
  expect_named(coef(cp), c("a1", "b1", "c1", "a2", "b2", "k"))
})

test_that("settings are enough to refit", {
  fit <- rc_powerlaw(discharge, stage, data = thompson, variance = "prop")
  refit <- do.call(rc_powerlaw, c(list(thompson$discharge, thompson$stage), fit$settings))
  expect_equal(coef(refit), coef(fit))
  fit2 <- rc_poly(discharge, stage, data = thompson, degree = 3)
  refit2 <- do.call(rc_poly, c(list(thompson$discharge, thompson$stage), fit2$settings))
  expect_equal(coef(refit2), coef(fit2))
})

test_that("the two-segment fit no longer takes conflev or predlev", {
  expect_error(
    rc_2seg_powerlaw(discharge, stage, data = thompson, kstart = 2, variance = "prop",
                conflev = 0.9),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("loess records its residual scale and equivalent parameters", {
  fit <- rc_loess(discharge, stage, data = thompson)
  expect_gt(fit$rse, 0)
  expect_gt(fit$enp, 1)
})

test_that("the fit records its weighting scheme, with any estimate", {
  fit <- rc_powerlaw(discharge, stage, data = thompson, variance = "power")
  expect_s3_class(fit$model, "gnls")
  expect_true(is.numeric(fit$variance$exponent) && length(fit$variance$exponent) == 1L)
  # settings keep the scheme as given, for refitting
  expect_null(fit$settings$variance$exponent)
  expect_s3_class(rc_powerlaw(discharge, stage, data = thompson)$variance, "rc_var_none")
  # the power describes the scatter, not the curve
  expect_named(coef(fit), c("a", "b", "c"))
  expect_named(fit$curve_parameters, c("a", "b", "c"))
})

# -- reweighting ---------------------------------------------------------------

test_that("converged reweighting is recorded", {
  for (fit in list(
    rc_powerlaw(discharge, stage, data = thompson, variance = "prop"),
    rc_poly(discharge, stage, data = thompson, variance = "prop"),
    rc_2seg_powerlaw(discharge, stage, data = thompson, variance = "prop", kstart = 2)
  )) {
    expect_true(fit$irls$converged)
    expect_gt(fit$irls$iterations, 1L)
    # the weights are those the final model was fitted with
    expect_equal(unname(fit$weights_used), unname(stats::weights(fit$model)))
  }
  expect_null(rc_powerlaw(discharge, stage, data = thompson)$irls)
})

test_that("reweighting that runs out of rounds warns and says so", {
  expect_warning(
    fit <- rc_powerlaw(discharge, stage, data = thompson, variance = var_prop(maxiter = 1)),
    "did not converge in 1 rounds"
  )
  expect_false(fit$irls$converged)
  expect_equal(fit$irls$iterations, 1L)
  expect_warning(
    rc_poly(discharge, stage, data = thompson, variance = var_prop(maxiter = 1)),
    "did not converge"
  )
  expect_warning(
    rc_2seg_powerlaw(discharge, stage, data = thompson, variance = var_prop(maxiter = 1), kstart = 2),
    "did not converge"
  )
})

test_that("reweighting stops on the change in fitted discharge", {
  fit <- rc_powerlaw(discharge, stage, data = thompson, variance = var_prop(tol = 1e-3))
  tight <- rc_powerlaw(discharge, stage, data = thompson, variance = var_prop(tol = 1e-10))
  expect_lte(fit$irls$iterations, tight$irls$iterations)
  expect_equal(stats::fitted(fit$model), stats::fitted(tight$model),
               tolerance = 1e-2)
})

test_that("every table the package returns is a tibble", {
  expect_s3_class(thompson, "tbl_df")
  for (nm in names(fits <- all_fits())) {
    expect_s3_class(fits[[nm]]$gaugings, "tbl_df")
    p <- suppressMessages(predict(fits[[nm]], new_stage = 3, conflev = 0.95))
    expect_s3_class(p, "tbl_df")
  }
  two <- fits$rc_2seg_powerlaw
  expect_s3_class(two$kstart_search, "tbl_df")
  expect_s3_class(
    suppressWarnings(boot_limits_2seg(two, new_stage = 3, conflev = 0.9, B = 5,
                                      seed = 1)),
    "tbl_df"
  )
})

test_that("loading the package loads tibble", {
  # without an import, tibble would load only on first use, and until then
  # the package's tibbles would behave like plain data frames
  expect_true("tibble" %in% names(getNamespaceImports("CSHShydRometry")))
})

test_that("no code relies on partial matching", {
  withr::local_options(
    warnPartialMatchDollar = TRUE,
    warnPartialMatchArgs = TRUE,
    warnPartialMatchAttr = TRUE
  )
  expect_no_warning({
    fits <- c(
      all_fits(),
      list(rc_powerlaw_log(discharge, stage, data = thompson,
                        offset = -1.3),
           rc_loess(discharge, stage, data = thompson, variance = "prop"))
    )
    for (fit in fits) {
      suppressMessages(
        predict(fit, new_stage = c(1, 3), conflev = 0.9, predlev = 0.9)
      )
      coef(fit)
    }
    boot_limits_2seg(fits$rc_2seg_powerlaw, new_stage = 2, conflev = 0.9,
                     predlev = 0.9, B = 3, seed = 1)
  })
})


test_that("fitted() and residuals() work for every fit", {
  fits <- all_fits()
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    n <- nrow(fit$gaugings)
    expect_length(fitted(fit), n)
    expect_length(residuals(fit), n)
    expect_true(all(is.finite(residuals(fit, divide_by_sd = TRUE))))
  }
  # on the scale the model is fitted on
  one <- fits$rc_powerlaw
  expect_equal(residuals(one), one$gaugings$discharge - fitted(one))
  log_fit <- fits$rc_powerlaw_log
  expect_equal(residuals(log_fit),
               log(log_fit$gaugings$discharge) - log(fitted(log_fit)))
  expect_error(residuals(one, divide_by_sd = "yes"))
  expect_error(fitted(fits[[1]], stage = 3), class = "rlib_error_dots_nonempty")
})

test_that("residuals divided by the sd match the underlying model's", {
  for (variance in list("none", "prop", "power")) {
    fit <- rc_powerlaw(discharge, stage, data = thompson, variance = variance)
    expect_equal(unname(fitted(fit)), unname(as.numeric(fitted(fit$model))))
    expect_equal(
      unname(residuals(fit, divide_by_sd = TRUE)),
      unname(as.numeric(residuals(fit$model, type = "pearson"))),
      tolerance = 1e-6
    )
  }
  log_fit <- rc_powerlaw_log(discharge, stage, data = thompson)
  expect_equal(
    unname(residuals(log_fit, divide_by_sd = TRUE)),
    as.numeric(residuals(log_fit$model)) / log_fit$rse
  )
})
