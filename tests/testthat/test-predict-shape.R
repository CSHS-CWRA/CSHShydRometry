# The promise this package makes about predict(): the columns you get back
# depend only on which of conflev and predlev you asked for -- never on the
# model, the weighting, or the interval method. That is what lets results from
# different approaches be stacked with rbind().

one_seg_fits <- function() {
  q <- thompson$q
  h <- thompson$h
  list(
    rc_log_ols = rc_log_ols(q, h),
    rc_log_nls = rc_log_nls(q, h),
    rc_nls = rc_nls(q, h),
    rc_gnls = rc_gnls(q, h),
    rc_poly = rc_poly(q, h),
    rc_loess = rc_loess(q, h)
  )
}

two_seg_fit <- function() {
  rc_nls_2seg(
    thompson$q,
    thompson$h,
    wts_code = "prop",
    kstart = 2
  )
}

test_that("every one-segment model returns the same columns", {
  fits <- one_seg_fits()
  want <- c("h", "fit", "ci_lwr", "ci_upr", "pi_lwr", "pi_upr")
  for (nm in names(fits)) {
    p <- suppressMessages(
      predict(fits[[nm]], hpred = c(1, 3), conflev = 0.95, predlev = 0.95)
    )
    expect_named(p, want, info = nm)
  }
})

test_that("results from different models stack with rbind()", {
  fits <- one_seg_fits()
  rows <- lapply(names(fits), function(nm) {
    p <- suppressMessages(
      predict(fits[[nm]], hpred = c(1, 3), conflev = 0.95, predlev = 0.95)
    )
    cbind(model = nm, as.data.frame(p))
  })
  stacked <- do.call(rbind, rows)
  expect_equal(nrow(stacked), 2 * length(fits))
})

test_that("asking for neither level returns neither set of columns", {
  fit <- rc_nls(thompson$q, thompson$h)
  p <- predict(fit, hpred = c(1, 3))
  expect_named(p, c("h", "fit"))
})

test_that("conflev alone returns only confidence columns", {
  fit <- rc_nls(thompson$q, thompson$h)
  p <- predict(fit, hpred = c(1, 3), conflev = 0.95)
  expect_named(p, c("h", "fit", "ci_lwr", "ci_upr"))
})

test_that("a single prediction stage works", {
  # regression test: indexing a one-row matrix used to drop it to a vector,
  # so predicting at exactly one stage failed
  fits <- one_seg_fits()
  for (nm in names(fits)) {
    p <- suppressMessages(
      predict(fits[[nm]], hpred = 3, conflev = 0.95, predlev = 0.95)
    )
    expect_equal(nrow(p), 1L, info = nm)
    expect_true(is.finite(p$fit), info = nm)
  }
  p2 <- predict(two_seg_fit(), hpred = 3, conflev = 0.95)
  expect_equal(nrow(p2), 1L)
})

test_that("specified weights give NA prediction limits, not missing columns", {
  keep <- !is.na(thompson$uq)
  d <- thompson[keep, ]
  skip_if(nrow(d) < 10, "too few gaugings with a reported uncertainty")
  fit <- rc_nls(d$q, d$h, wts_code = "spec", wts = 1 / d$uq^2)
  p <- suppressMessages(
    predict(fit, hpred = c(1, 3), conflev = 0.95, predlev = 0.95)
  )
  expect_true(all(c("pi_lwr", "pi_upr") %in% names(p)))
  expect_true(all(is.na(p$pi_lwr)))
  expect_true(all(is.finite(p$ci_lwr)))
})
