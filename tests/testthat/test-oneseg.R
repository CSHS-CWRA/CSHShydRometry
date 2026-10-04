# The one-segment constructors: input handling, weighting schemes, and the
# limits each weighting leads to.

# The gaugings with a reported uncertainty, and variances from it: the
# uncertainty is a percentage of the discharge at two standard deviations.
# The gauging reported as 0.0226% is left out: it is probably a fraction
# entered as a percentage (see ?thompson), and its weight would be over 99.9%
# of the total, forcing the curve through it.
spec_data <- function() {
  d <- thompson[!is.na(thompson$uncertainty_pct) & thompson$uncertainty_pct > 1, ]
  d$variance <- (d$uncertainty_pct / 100 * d$discharge / 2)^2
  skip_if(nrow(d) < 10, "too few gaugings with a reported uncertainty")
  d
}

test_that("discharge and stage may be vectors or columns of data", {
  from_vectors <- rc_powerlaw(thompson$discharge, thompson$stage)
  from_columns <- rc_powerlaw(discharge, stage, data = thompson)
  expect_equal(from_vectors$pars, from_columns$pars)
})

test_that("arguments after discharge and stage must be named", {
  expect_error(rc_powerlaw(discharge, stage, thompson), class = "rlib_error_dots_nonempty")
  expect_error(
    rc_poly(discharge, stage, data = thompson, degre = 3),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("column names are found when called from inside another function", {
  wrapper <- function(d) rc_powerlaw_log(discharge, stage, data = d)
  expect_equal(wrapper(thompson)$pars, rc_powerlaw_log(discharge, stage, data = thompson)$pars)
  local_q <- thompson$discharge
  local_h <- thompson$stage
  inner <- function() rc_powerlaw_log(local_q, local_h)
  expect_s3_class(inner(), "rc_powerlaw_log")
})

test_that("gaugings with a missing stage or discharge are dropped", {
  d <- thompson
  d$discharge[c(3, 10)] <- NA
  d$stage[20] <- NA
  fit <- rc_powerlaw(discharge, stage, data = d)
  expect_equal(nrow(fit$gaugings), nrow(thompson) - 3L)
  expect_false(anyNA(fit$gaugings))
})

test_that("var_spec() takes a vector, not a column of the fit's data", {
  d <- spec_data()
  d$uncertainty_sd <- d$uncertainty_pct / 100 * d$discharge / 2
  # evaluated where it is called, so a bare column name is not found
  expect_error(var_spec(uncertainty_sd^2), "uncertainty_sd")
  expect_error(var_spec("not numbers"))
  fits <- list(
    rc_powerlaw(discharge, stage, data = d, variance = var_spec(d$uncertainty_sd^2)),
    rc_poly(discharge, stage, data = d, variance = var_spec(d$uncertainty_sd^2)),
    rc_loess(discharge, stage, data = d, variance = var_spec(d$uncertainty_sd^2),
             span = 1)
  )
  for (fit in fits) {
    expect_equal(fit$settings$variance$values, d$uncertainty_sd^2)
  }
})

test_that("specified variances stay aligned when gaugings are dropped", {
  d <- spec_data()
  w <- d$variance
  d$discharge[2] <- NA
  fit <- rc_powerlaw(discharge, stage, data = d, variance = var_spec(w))
  expect_equal(fit$weights_used, 1 / w[-2])
})

test_that("weights of the wrong length are an error", {
  expect_error(rc_powerlaw(discharge, stage, data = thompson, variance = var_spec(1:3)))
  expect_error(rc_poly(discharge, stage, data = thompson, variance = var_spec(1:3)))
  expect_error(rc_loess(discharge, stage, data = thompson, variance = var_spec(1:3)))
})

test_that("invalid levels are rejected", {
  fit <- rc_powerlaw(discharge, stage, data = thompson)
  expect_error(predict(fit, new_stage = 3, conflev = 1.5))
  expect_error(predict(fit, new_stage = 3, predlev = -0.1))
})

test_that("an unknown weighting scheme is an error", {
  expect_error(rc_powerlaw(discharge, stage, data = thompson, variance = "bogus"))
})

test_that("predict() defaults to the observed stage range", {
  fit <- rc_powerlaw_log(discharge, stage, data = thompson)
  p <- predict(fit)
  expect_equal(nrow(p), 1000L)
  expect_equal(range(p$stage), range(thompson$stage))
})

test_that("print() describes the fit and returns it invisibly", {
  fit <- rc_powerlaw(discharge, stage, data = thompson)
  expect_output(out <- withVisible(print(fit)), "Method: rc_powerlaw")
  expect_false(out$visible)
  expect_identical(out$value, fit)
})

# -- rc_powerlaw --------------------------------------------------------------------

test_that("rc_powerlaw: proportional weights give sensible limits", {
  fit <- rc_powerlaw(discharge, stage, data = thompson, variance = "prop")
  expect_equal(fit$settings$variance$type, "prop")
  p <- predict(fit, new_stage = c(1, 3, 6), conflev = 0.95, predlev = 0.95)
  expect_true(all(p$ci_lwr < p$fit & p$fit < p$ci_upr))
  expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr))
  # constant coefficient of variation: the prediction band widens with flow
  expect_true(all(diff(p$pi_upr - p$pi_lwr) > 0))
})

test_that("rc_powerlaw: specified weights give NA prediction limits with a note", {
  d <- spec_data()
  fit <- rc_powerlaw(discharge, stage, data = d, variance = var_spec(d$variance))
  expect_message(
    p <- predict(fit, new_stage = c(1, 3), conflev = 0.95, predlev = 0.95),
    "specified weights"
  )
  expect_true(all(is.finite(p$ci_lwr)))
  expect_true(all(is.na(p$pi_upr)))
})

# -- rc_poly -------------------------------------------------------------------

test_that("rc_poly: degree sets the number of coefficients", {
  for (deg in 1:3) {
    fit <- rc_poly(discharge, stage, data = thompson, degree = deg)
    expect_named(fit$pars, paste0("b", 0:deg))
  }
  expect_error(rc_poly(discharge, stage, data = thompson, degree = 0))
})

test_that("rc_poly: a straight line matches lm()", {
  fit <- rc_poly(discharge, stage, data = thompson, degree = 1)
  ref <- stats::lm(discharge ~ stage, data = thompson)
  expect_equal(unname(coef(fit)), unname(stats::coef(ref)))
})

test_that("rc_poly: proportional weights fit and give limits", {
  fit <- rc_poly(discharge, stage, data = thompson, variance = "prop")
  expect_s3_class(fit$model, "nls")
  p <- predict(fit, new_stage = c(1, 3, 6), conflev = 0.95, predlev = 0.95)
  expect_named(p, c("stage", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr"))
  expect_true(all(p$ci_lwr < p$fit & p$fit < p$ci_upr))
  expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr))
})

test_that("rc_poly: specified weights give NA prediction limits with a note", {
  d <- spec_data()
  fit <- rc_poly(discharge, stage, data = d, variance = var_spec(d$variance))
  expect_equal(fit$settings$variance$type, "spec")
  expect_message(
    p <- predict(fit, new_stage = c(1, 3), conflev = 0.95, predlev = 0.95),
    "specified weights"
  )
  expect_true(all(is.finite(p$ci_lwr)))
  expect_true(all(is.na(p$pi_lwr)))
})

# -- rc_loess ------------------------------------------------------------------

test_that("rc_loess: every weighting fits", {
  d <- spec_data()
  fits <- list(
    none = rc_loess(discharge, stage, data = thompson),
    prop = rc_loess(discharge, stage, data = thompson, variance = "prop"),
    spec = rc_loess(discharge, stage, data = d, variance = var_spec(d$variance),
                    span = 1)
  )
  for (nm in names(fits)) {
    p <- predict(fits[[nm]], new_stage = 3, conflev = 0.95)
    expect_true(is.finite(p$fit), info = nm)
    expect_true(p$ci_lwr < p$ci_upr, info = nm)
  }
})

test_that("rc_loess: without extrapolation, stages outside the data give NA", {
  fit <- rc_loess(discharge, stage, data = thompson, extrapolate = FALSE)
  p <- predict(fit, new_stage = c(max(thompson$stage) + 1, 3))
  expect_true(is.na(p$fit[1]))
  expect_true(is.finite(p$fit[2]))
})

test_that("rc_loess: prediction limits are NA with a note", {
  fit <- rc_loess(discharge, stage, data = thompson)
  expect_message(
    p <- predict(fit, new_stage = 3, predlev = 0.95),
    "not implemented"
  )
  expect_true(is.na(p$pi_lwr))
})

# -- log scale -----------------------------------------------------------------

test_that("log-scale fits report bias-corrected coefficients", {
  for (fit in list(rc_powerlaw_log(discharge, stage, data = thompson),
                   rc_powerlaw_log(discharge, stage, data = thompson, offset = -1.3))) {
    expect_gt(fit$a_corrected[["nbc"]], fit$pars$a)
    expect_true(fit$pars$c < min(thompson$stage))
  }
})

test_that("log-scale limits are positive", {
  for (fit in list(rc_powerlaw_log(discharge, stage, data = thompson),
                   rc_powerlaw_log(discharge, stage, data = thompson, offset = -1.3))) {
    p <- predict(fit, new_stage = c(1, 3), conflev = 0.95, predlev = 0.95)
    expect_true(all(p$pi_lwr > 0))
  }
})

test_that("a known c is held fixed, and gives a linear fit on the log scale", {
  fit <- rc_powerlaw_log(discharge, stage, data = thompson, offset = -1.3)
  expect_equal(fit$pars$c, -1.3)
  expect_equal(fit$settings$offset, -1.3)
  expect_s3_class(fit$model, "lm")
  expect_null(rc_powerlaw_log(discharge, stage, data = thompson)$settings$offset)
  expect_error(
    rc_powerlaw_log(discharge, stage, data = thompson, offset = min(thompson$stage)),
    "`offset` must be below"
  )
})

test_that("fixing c at its estimate keeps the curve but narrows the limits", {
  # this is what the former rc_log_ols() did: estimate c, then treat it as
  # known, leaving its uncertainty out of the limits
  est <- rc_powerlaw_log(discharge, stage, data = thompson)
  fixed <- rc_powerlaw_log(discharge, stage, data = thompson, offset = est$pars$c)
  expect_equal(fixed$pars$a, est$pars$a, tolerance = 1e-6)
  expect_equal(fixed$pars$b, est$pars$b, tolerance = 1e-6)
  p_est <- predict(est, new_stage = 3, conflev = 0.95)
  p_fixed <- predict(fixed, new_stage = 3, conflev = 0.95)
  expect_equal(p_fixed$fit, p_est$fit, tolerance = 1e-6)
  expect_lt(p_fixed$ci_upr - p_fixed$ci_lwr, p_est$ci_upr - p_est$ci_lwr)
})

test_that("the log-scale fit explains why it takes no variance scheme", {
  expect_error(
    rc_powerlaw_log(discharge, stage, data = thompson, variance = "prop"),
    "not implemented for the log-scale fit"
  )
  # other stray arguments still get the usual error
  expect_error(
    rc_powerlaw_log(discharge, stage, data = thompson, zero_flow = 1),
    class = "rlib_error_dots_nonempty"
  )
})


test_that("rc_powerlaw() holds a given offset fixed, under every weighting", {
  for (scheme in list("none", "prop", "power")) {
    est <- rc_powerlaw(discharge, stage, data = thompson, variance = scheme)
    fixed <- rc_powerlaw(discharge, stage, data = thompson, variance = scheme,
                         offset = est$pars$c)
    expect_equal(fixed$pars$c, est$pars$c)
    expect_equal(fixed$settings$offset, est$pars$c)
    expect_named(coef(fixed$model), c("a", "b"))
    # the same curve, with limits that leave out the uncertainty in c
    expect_equal(coef(fixed), coef(est), tolerance = 1e-4)
    p_est <- predict(est, new_stage = 3, conflev = 0.95, predlev = 0.95)
    p_fixed <- predict(fixed, new_stage = 3, conflev = 0.95, predlev = 0.95)
    expect_equal(p_fixed$fit, p_est$fit, tolerance = 1e-4)
    expect_lt(p_fixed$ci_upr - p_fixed$ci_lwr, p_est$ci_upr - p_est$ci_lwr)
  }
  expect_error(
    rc_powerlaw(discharge, stage, data = thompson, offset = min(thompson$stage)),
    "`offset` must be below"
  )
})
