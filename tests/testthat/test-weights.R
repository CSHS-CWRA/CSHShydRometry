# Weighting schemes, and the settings passed to nls().

test_that("the shorthand strings match the constructors", {
  a <- rc_powerlaw(discharge, stage, data = thompson, wts = "prop")
  b <- rc_powerlaw(discharge, stage, data = thompson, wts = wts_prop())
  expect_equal(coef(a), coef(b))
  expect_equal(
    coef(rc_powerlaw(discharge, stage, data = thompson, wts = "none")),
    coef(rc_powerlaw(discharge, stage, data = thompson))
  )
})

test_that("invalid weighting schemes are errors", {
  expect_error(
    rc_powerlaw(discharge, stage, data = thompson, wts = "spec"),
    "wts_spec"
  )
  expect_error(rc_powerlaw(discharge, stage, data = thompson, wts = "bogus"))
  expect_error(rc_powerlaw(discharge, stage, data = thompson, wts = 1:3), "must be")
  expect_error(wts_prop(tol = -1))
  expect_error(wts_prop(maxiter = 0))
})

test_that("the fit records its weighting scheme", {
  fit <- rc_powerlaw(discharge, stage, data = thompson, wts = wts_prop(tol = 1e-4))
  expect_s3_class(fit$settings$wts, "rc_wts_prop")
  expect_equal(fit$settings$wts$tol, 1e-4)
  expect_equal(fit$settings$wts$maxiter, 100)
})

test_that("weighting schemes print what they are", {
  expect_output(print(wts_none()), "none")
  expect_output(print(wts_prop(maxiter = 20)), "maxiter = 20")
  expect_output(print(wts_spec(c(1, 2, 3))), "for 3 gaugings")
  # leave out the doubtful 0.0226% gauging; see ?thompson
  d <- thompson[!is.na(thompson$uncertainty_pct) & thompson$uncertainty_pct > 1, ]
  d$uncertainty_sd <- d$uncertainty_pct / 100 * d$discharge / 2
  fit <- rc_powerlaw(discharge, stage, data = d,
                wts = wts_spec(1 / d$uncertainty_sd^2))
  expect_output(print(fit$settings$wts), "for 18 gaugings")
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
    rc_powerlaw(discharge, stage, data = thompson, control = tight),
    "maximum of 1"
  )
  expect_error(
    rc_powerlaw_log(discharge, stage, data = thompson, control = tight),
    "maximum of 1"
  )
  expect_error(
    rc_2seg_powerlaw(discharge, stage, data = thompson, wts = "prop",
                control = list(maxiter = 1)),
    "failed from every starting breakpoint"
  )
  expect_error(
    rc_powerlaw(discharge, stage, data = thompson, control = list(1)),
    "names"
  )
})

test_that("the old arguments are gone", {
  for (old in list(
    list(wts_code = "prop"), list(wts_tol = 1e-3), list(wts_maxiter = 5),
    list(nls_tol = 1e-3), list(nls_maxiter = 5), list(config = "add")
  )) {
    expect_error(
      do.call(rc_2seg_powerlaw, c(list(thompson$discharge, thompson$stage), old)),
      class = "rlib_error_dots_nonempty"
    )
  }
})

# -- wts_power(): the power of the flow estimated -------------------------------

test_that("wts_power() estimates the power with nlme::gnls()", {
  fit <- rc_powerlaw(discharge, stage, data = thompson, wts = wts_power())
  # the same as nlme::gnls() from the same start
  d <- as.data.frame(thompson[, c("discharge", "stage")])
  ref <- nlme::gnls(
    discharge ~ a * (stage - c)^b,
    data = d,
    start = as.list(coef(rc_powerlaw(discharge, stage, data = thompson,
                                  wts = "prop"))),
    weights = nlme::varPower(),
    control = nlme::gnlsControl(maxIter = 1e5, minScale = 1e-5)
  )
  expect_equal(coef(fit), coef(ref))
  expect_equal(
    fit$wts$exponent,
    unname(coef(ref$modelStruct$varStruct, unconstrained = FALSE))
  )
  expect_equal(coef(rc_powerlaw(discharge, stage, data = thompson, wts = "power")),
               coef(fit))
  p <- predict(fit, new_stage = c(1, 3), conflev = 0.95, predlev = 0.95)
  expect_named(p, c("stage", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr"))
  expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr))
})

test_that("wts_power() is for rc_powerlaw() only, and says why", {
  for (fitter in list(rc_poly, rc_loess)) {
    expect_error(
      fitter(discharge, stage, data = thompson, wts = wts_power()),
      "available in `rc_powerlaw\\(\\)` only.*nlme::gnls"
    )
  }
  expect_error(
    rc_2seg_powerlaw(discharge, stage, data = thompson, wts = "power"),
    "not in `rc_2seg_powerlaw\\(\\)`"
  )
})

test_that("the schemes print what they are", {
  expect_output(print(wts_prop()), "proportional to the flow \\(tol")
  expect_output(print(wts_power()), "power to be estimated")
  fit <- rc_powerlaw(discharge, stage, data = thompson, wts = wts_power())
  expect_output(print(fit$wts), "to the power [0-9.]+ \\(estimated\\)")
})

test_that("both schemes record nlme's description of the scatter", {
  expect_s3_class(wts_prop()$variance, "varPower")
  expect_s3_class(wts_power()$variance, "varPower")
})

test_that("wts_power() limits are reproducible", {
  fit <- rc_powerlaw(discharge, stage, data = thompson, wts = wts_power())
  at <- c(1, 4, 8)
  expect_identical(
    predict(fit, new_stage = at, conflev = 0.95, predlev = 0.95),
    predict(fit, new_stage = at, conflev = 0.95, predlev = 0.95)
  )
})

test_that("wts_power() limits at a power of 1 match those under wts_prop()", {
  # gnls() with the power fixed at 1 reaches the same fit as the rounds
  prop <- rc_powerlaw(discharge, stage, data = thompson, wts = wts_prop(tol = 1e-10))
  ref <- nlme::gnls(
    discharge ~ a * (stage - c)^b,
    data = as.data.frame(thompson[, c("discharge", "stage")]),
    start = as.list(coef(prop)),
    weights = nlme::varPower(fixed = 1)
  )
  at <- c(1, 4, 8)
  expect_equal(
    gnls_limits(ref, data.frame(stage = at), pars = as.list(coef(ref)), power = 1, conflev = 0.9,
                predlev = 0.9),
    predict(prop, new_stage = at, conflev = 0.9, predlev = 0.9),
    tolerance = 1e-5
  )
})
