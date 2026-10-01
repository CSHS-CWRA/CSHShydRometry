# The promise this package makes about predict(): the columns you get back
# depend only on which of conflev and predlev you asked for -- never on the
# model, the weighting, or the interval method. That is what lets results from
# different approaches be stacked with rbind().

one_seg_fits <- function() {
  discharge <- thompson$discharge
  stage <- thompson$stage
  list(
    rc_power_log_c = rc_power_log(discharge, stage, zero_flow_stage = -1.3),
    rc_power_log = rc_power_log(discharge, stage),
    rc_power = rc_power(discharge, stage),
    rc_power_wts_power = rc_power(discharge, stage, wts = wts_power()),
    rc_poly = rc_poly(discharge, stage),
    rc_loess = rc_loess(discharge, stage)
  )
}

two_seg_fit <- function() {
  rc_2seg_power(
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
      predict(fits[[nm]], stage = c(1, 3), conflev = 0.95, predlev = 0.95)
    )
    expect_named(p, want, info = nm)
  }
})

test_that("results from different models stack with rbind()", {
  fits <- one_seg_fits()
  rows <- lapply(names(fits), function(nm) {
    p <- suppressMessages(
      predict(fits[[nm]], stage = c(1, 3), conflev = 0.95, predlev = 0.95)
    )
    cbind(model = nm, as.data.frame(p))
  })
  stacked <- do.call(rbind, rows)
  expect_equal(nrow(stacked), 2 * length(fits))
})

test_that("asking for neither level returns neither set of columns", {
  fit <- rc_power(thompson$discharge, thompson$stage)
  p <- predict(fit, stage = c(1, 3))
  expect_named(p, c("stage", "fit"))
})

test_that("conflev alone returns only confidence columns", {
  fit <- rc_power(thompson$discharge, thompson$stage)
  p <- predict(fit, stage = c(1, 3), conflev = 0.95)
  expect_named(p, c("stage", "fit", "ci_lwr", "ci_upr"))
})

test_that("a single prediction stage works", {
  # regression test: indexing a one-row matrix used to drop it to a vector,
  # so predicting at exactly one stage failed
  fits <- one_seg_fits()
  for (nm in names(fits)) {
    p <- suppressMessages(
      predict(fits[[nm]], stage = 3, conflev = 0.95, predlev = 0.95)
    )
    expect_equal(nrow(p), 1L, info = nm)
    expect_true(is.finite(p$fit), info = nm)
  }
  p2 <- predict(two_seg_fit(), stage = 3, conflev = 0.95)
  expect_equal(nrow(p2), 1L)
})

test_that("specified weights give NA prediction limits, not missing columns", {
  keep <- !is.na(thompson$uncertainty_pct)
  d <- thompson[keep, ]
  skip_if(nrow(d) < 10, "too few gaugings with a reported uncertainty")
  sd <- d$uncertainty_pct / 100 * d$discharge / 2
  fit <- rc_power(d$discharge, d$stage, wts = wts_spec(1 / sd^2))
  p <- suppressMessages(
    predict(fit, stage = c(1, 3), conflev = 0.95, predlev = 0.95)
  )
  expect_true(all(c("pi_lwr", "pi_upr") %in% names(p)))
  expect_true(all(is.na(p$pi_lwr)))
  expect_true(all(is.finite(p$ci_lwr)))
})
