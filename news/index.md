# Changelog

## CSHShydRometry (development version)

- First CRAN release.

- Fitting functions are named for the model they fit, such as
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md),
  [`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md)
  and
  [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md).
  Their arguments are `discharge` and `stage` (formerly `q` and `h`),
  and every argument after those must be named.

- How the gaugings scatter about the curve is chosen by a single
  `variance` argument taking a variance scheme.

- `c`, the offset, can be held at a known value with `offset` in
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md)
  and
  [`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md).
  `rc_log_ols()` is removed.

- [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md)
  uses `combine = "replace"` or `"add"` to say how the segments combine
  above the breakpoint (formerly `config`, `"piecewise"` or
  `"compound"`), and tries several starting breakpoints and keeps the
  best fit.

- Bug fixes: `rc_poly(degree = 1)` fitted spurious higher-order terms,
  and user-supplied weights fell out of step with the gaugings when any
  had a missing value.

## CSHShydRometry 0.0.1

- First packaged version, a prototype, not released. Derived from
  rating-curve functions written by R. Dan Moore.
