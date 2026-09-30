# The one-segment constructors: input handling, weighting schemes, and the
# limits each weighting leads to.

# The gaugings with a reported uncertainty, and weights from it: the
# uncertainty is a percentage of the discharge at two standard deviations.
spec_data <- function() {
  d <- thompson[!is.na(thompson$uncertainty_pct), ]
  d$wts <- 1 / (d$uncertainty_pct / 100 * d$discharge / 2)^2
  skip_if(nrow(d) < 10, "too few gaugings with a reported uncertainty")
  d
}

test_that("discharge and stage may be vectors or columns of data", {
  from_vectors <- rc_nls(thompson$discharge, thompson$stage)
  from_columns <- rc_nls(discharge, stage, data = thompson)
  expect_equal(from_vectors$pars, from_columns$pars)
})

test_that("arguments after discharge and stage must be named", {
  expect_error(rc_nls(discharge, stage, thompson), class = "rlib_error_dots_nonempty")
  expect_error(
    rc_poly(discharge, stage, data = thompson, degre = 3),
    class = "rlib_error_dots_nonempty"
  )
})

test_that("column names are found when called from inside another function", {
  wrapper <- function(d) rc_log_ols(discharge, stage, data = d)
  expect_equal(wrapper(thompson)$pars, rc_log_ols(discharge, stage, data = thompson)$pars)
  local_q <- thompson$discharge
  local_h <- thompson$stage
  inner <- function() rc_log_nls(local_q, local_h)
  expect_s3_class(inner(), "rc_log_nls")
})

test_that("gaugings with a missing stage or discharge are dropped", {
  d <- thompson
  d$discharge[c(3, 10)] <- NA
  d$stage[20] <- NA
  fit <- rc_gnls(discharge, stage, data = d)
  expect_equal(nrow(fit$gaugings), nrow(thompson) - 3L)
  expect_false(anyNA(fit$gaugings))
})

test_that("specified weights stay aligned when gaugings are dropped", {
  d <- spec_data()
  w <- d$wts
  d$discharge[2] <- NA
  fit <- rc_nls(discharge, stage, data = d, wts_code = "spec", wts = w)
  expect_equal(fit$weights, w[-2])
})

test_that("weights of the wrong length are an error", {
  expect_error(rc_nls(discharge, stage, data = thompson, wts_code = "spec", wts = 1:3))
  expect_error(rc_poly(discharge, stage, data = thompson, wts_code = "spec", wts = 1:3))
  expect_error(rc_loess(discharge, stage, data = thompson, wts_code = "spec", wts = 1:3))
})

test_that("invalid levels are rejected", {
  fit <- rc_nls(discharge, stage, data = thompson)
  expect_error(predict(fit, stage = 3, conflev = 1.5))
  expect_error(predict(fit, stage = 3, predlev = -0.1))
})

test_that("an unknown weighting scheme is an error", {
  expect_error(rc_nls(discharge, stage, data = thompson, wts_code = "bogus"))
})

test_that("predict() defaults to the observed stage range", {
  fit <- rc_log_ols(discharge, stage, data = thompson)
  p <- predict(fit)
  expect_equal(nrow(p), 1000L)
  expect_equal(range(p$stage), range(thompson$stage))
})

test_that("print() describes the fit and returns it invisibly", {
  fit <- rc_nls(discharge, stage, data = thompson)
  expect_output(out <- withVisible(print(fit)), "Method: rc_nls")
  expect_false(out$visible)
  expect_identical(out$value, fit)
})

# -- rc_nls --------------------------------------------------------------------

test_that("rc_nls: proportional weights give sensible limits", {
  fit <- rc_nls(discharge, stage, data = thompson, wts_code = "prop")
  expect_equal(fit$settings$wts_code, "prop")
  p <- predict(fit, stage = c(1, 3, 6), conflev = 0.95, predlev = 0.95)
  expect_true(all(p$ci_lwr < p$fit & p$fit < p$ci_upr))
  expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr))
  # constant coefficient of variation: the prediction band widens with flow
  expect_true(all(diff(p$pi_upr - p$pi_lwr) > 0))
})

test_that("rc_nls: specified weights give NA prediction limits with a note", {
  d <- spec_data()
  fit <- rc_nls(discharge, stage, data = d, wts_code = "spec", wts = d$wts)
  expect_message(
    p <- predict(fit, stage = c(1, 3), conflev = 0.95, predlev = 0.95),
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
  fit <- rc_poly(discharge, stage, data = thompson, wts_code = "prop")
  expect_s3_class(fit$model, "nls")
  p <- predict(fit, stage = c(1, 3, 6), conflev = 0.95, predlev = 0.95)
  expect_named(p, c("stage", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr"))
  expect_true(all(p$ci_lwr < p$fit & p$fit < p$ci_upr))
  expect_true(all(p$pi_lwr < p$ci_lwr & p$ci_upr < p$pi_upr))
})

test_that("rc_poly: specified weights give NA prediction limits with a note", {
  d <- spec_data()
  fit <- rc_poly(discharge, stage, data = d, wts_code = "spec", wts = d$wts)
  expect_equal(fit$settings$wts_code, "spec")
  expect_message(
    p <- predict(fit, stage = c(1, 3), conflev = 0.95, predlev = 0.95),
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
    prop = rc_loess(discharge, stage, data = thompson, wts_code = "prop"),
    spec = rc_loess(discharge, stage, data = d, wts_code = "spec", wts = d$wts,
                    span = 1)
  )
  for (nm in names(fits)) {
    p <- predict(fits[[nm]], stage = 3, conflev = 0.95)
    expect_true(is.finite(p$fit), info = nm)
    expect_true(p$ci_lwr < p$ci_upr, info = nm)
  }
})

test_that("rc_loess: without extrapolation, stages outside the data give NA", {
  fit <- rc_loess(discharge, stage, data = thompson, extrapolate = FALSE)
  p <- predict(fit, stage = c(max(thompson$stage) + 1, 3))
  expect_true(is.na(p$fit[1]))
  expect_true(is.finite(p$fit[2]))
})

test_that("rc_loess: prediction limits are NA with a note", {
  fit <- rc_loess(discharge, stage, data = thompson)
  expect_message(
    p <- predict(fit, stage = 3, predlev = 0.95),
    "not implemented"
  )
  expect_true(is.na(p$pi_lwr))
})

# -- log scale -----------------------------------------------------------------

test_that("log-scale fits report bias-corrected coefficients", {
  for (fit in list(rc_log_ols(discharge, stage, data = thompson),
                   rc_log_nls(discharge, stage, data = thompson))) {
    expect_gt(fit$a_corrected[["nbc"]], fit$pars$a)
    expect_true(fit$pars$c < min(thompson$stage))
  }
})

test_that("log-scale limits are positive", {
  for (fit in list(rc_log_ols(discharge, stage, data = thompson),
                   rc_log_nls(discharge, stage, data = thompson))) {
    p <- predict(fit, stage = c(1, 3), conflev = 0.95, predlev = 0.95)
    expect_true(all(p$pi_lwr > 0))
  }
})
