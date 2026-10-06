# CSHShydRometry (development version)

Not yet released. The major changes since the first packaged version:

* Fitting functions are named for the model they fit: `rc_powerlaw()`,
  `rc_powerlaw_log()` and `rc_2seg_powerlaw()` (formerly `rc_nls()`,
  `rc_log_nls()` and `rc_nls_2seg()`), alongside `rc_poly()` and
  `rc_loess()`. Their arguments are `discharge` and `stage` (formerly `q` and
  `h`), and every argument after those must be named.

* How the gaugings scatter about the curve is chosen by a single `variance`
  argument taking a variance scheme: `var_none()`, `var_prop()`,
  `var_power()` or `var_spec()`; see `?variance`. It replaces `wts_code` and
  its settings. `rc_gnls()` is now `variance = var_power()` in
  `rc_powerlaw()`. This interface is experimental.

* `c`, the offset, can be held at a known value with `offset` in
  `rc_powerlaw()` and `rc_powerlaw_log()`. `rc_log_ols()` is removed.

* `rc_2seg_powerlaw()`:
  * `combine = "replace"` or `"add"` says how the segments combine above the
    breakpoint (formerly `config`, `"piecewise"` or `"compound"`);
  * it tries several starting breakpoints and keeps the best fit;
  * limits come from `predict(method = "delta")` or
    `predict(method = "boot")`. The `"sim"` method is removed, as its limits
    were often meaningless.

* `predict()` takes the stages as `new_stage`, returns a tibble with the
  same columns for every model, and gives an error for any argument it does
  not use. Limits under `var_power()` are computed by the delta method
  rather than simulated, and nlraa is no longer a dependency.

* New `fitted()` and `residuals()` methods. Fit objects are restructured
  the same way for every model; see `?rating_curve`.

* Two vignettes: `vignette("fitting")` and `vignette("uncertainty")`.

* Bug fixes: `rc_poly(degree = 1)` fitted spurious higher-order terms, and
  user-supplied weights fell out of step with the gaugings when any had a
  missing value.

# CSHShydRometry 0.0.1

* First packaged version, not released. Derived from rating-curve functions
  written by R. Dan Moore.
