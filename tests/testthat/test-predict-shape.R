# The promise this package makes about predict(): the columns you get back
# depend only on which of conflev and predlev you asked for -- never on the
# model, the weighting, or the interval method. That is what lets results from
# different approaches be stacked with rbind().

one_seg_fits <- function() {
  discharge <- thompson$discharge
  stage <- thompson$stage
  list(
    rc_powerlaw_log_offset = rc_powerlaw_log(discharge, stage, offset = -1.3),
    rc_powerlaw_log = rc_powerlaw_log(discharge, stage),
    rc_powerlaw = rc_powerlaw(discharge, stage),
    rc_power_wts_power = rc_powerlaw(discharge, stage, wts = wts_power()),
    rc_poly = rc_poly(discharge, stage),
    rc_loess = rc_loess(discharge, stage)
  )
}

two_seg_fit <- function() {
  rc_2seg_powerlaw(
    thompson$discharge,
    thompson$stage,
    wts = "prop",
    kstart = 2
  )
}

test_that("every one-segment model returns the same columns", {
  fits <- one_seg_fits()
  want <- c("stage", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr")
  for (nm in names(fits)) {
    p <- suppressMessages(
      predict(fits[[nm]], new_stage = c(1, 3), conflev = 0.95, predlev = 0.95)
    )
    expect_named(p, want, info = nm)
  }
})

test_that("results from different models stack with rbind()", {
  fits <- one_seg_fits()
  rows <- lapply(names(fits), function(nm) {
    p <- suppressMessages(
      predict(fits[[nm]], new_stage = c(1, 3), conflev = 0.95, predlev = 0.95)
    )
    cbind(model = nm, as.data.frame(p))
  })
  stacked <- do.call(rbind, rows)
  expect_equal(nrow(stacked), 2 * length(fits))
})

test_that("asking for neither level returns neither set of columns", {
  fit <- rc_powerlaw(thompson$discharge, thompson$stage)
  p <- predict(fit, new_stage = c(1, 3))
  expect_named(p, c("stage", "fit"))
})

test_that("conflev alone returns only confidence columns", {
  fit <- rc_powerlaw(thompson$discharge, thompson$stage)
  p <- predict(fit, new_stage = c(1, 3), conflev = 0.95)
  expect_named(p, c("stage", "fit", "ci_lwr", "ci_upr"))
})

test_that("a single prediction stage works", {
  # regression test: indexing a one-row matrix used to drop it to a vector,
  # so predicting at exactly one stage failed
  fits <- one_seg_fits()
  for (nm in names(fits)) {
    p <- suppressMessages(
      predict(fits[[nm]], new_stage = 3, conflev = 0.95, predlev = 0.95)
    )
    expect_equal(nrow(p), 1L, info = nm)
    expect_true(is.finite(p$fit), info = nm)
  }
  p2 <- predict(two_seg_fit(), new_stage = 3, conflev = 0.95)
  expect_equal(nrow(p2), 1L)
})

test_that("specified weights give NA prediction limits, not missing columns", {
  # leave out the doubtful 0.0226% gauging; see ?thompson
  keep <- !is.na(thompson$uncertainty_pct) & thompson$uncertainty_pct > 1
  d <- thompson[keep, ]
  skip_if(nrow(d) < 10, "too few gaugings with a reported uncertainty")
  sd <- d$uncertainty_pct / 100 * d$discharge / 2
  fit <- rc_powerlaw(d$discharge, d$stage, wts = wts_spec(1 / sd^2))
  p <- suppressMessages(
    predict(fit, new_stage = c(1, 3), conflev = 0.95, predlev = 0.95)
  )
  expect_true(all(c("pi_lwr", "pi_upr") %in% names(p)))
  expect_true(all(is.na(p$pi_lwr)))
  expect_true(all(is.finite(p$ci_lwr)))
})

test_that("predict() rejects arguments it does not take", {
  one <- rc_powerlaw(discharge, stage, data = thompson)
  for (fit in list(
    one,
    rc_powerlaw_log(discharge, stage, data = thompson),
    rc_poly(discharge, stage, data = thompson),
    rc_loess(discharge, stage, data = thompson)
  )) {
    expect_error(predict(fit, new_stage = 3, method = "boot"),
                 class = "rlib_error_dots_nonempty")
  }
  # the old name for new_stage, and a misspelt level
  expect_error(predict(one, stage = 3), class = "rlib_error_dots_nonempty")
  expect_error(predict(one, new_stage = 3, conflevel = 0.9),
               class = "rlib_error_dots_nonempty")
  # bootstrap settings with the delta method
  two <- two_seg_fit()
  expect_error(predict(two, new_stage = 3, B = 50),
               class = "rlib_error_dots_nonempty")
})
