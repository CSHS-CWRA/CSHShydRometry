# CSHShydRometry (development version)

## Breaking changes

* `rc_nls_2seg()` is renamed `rc_2seg_nls()`, and its class `rc_nls_2seg`
  is now `rc_2seg_nls`. Functions are named by number of segments first:
  one-segment curves keep their names (`rc_nls()`, `rc_poly()`, ...), and
  two-segment curves are `rc_2seg_*()`, leaving room for more two-segment
  methods and for curves with more segments.

* Every argument after the mandatory ones (`q` and `h`, or `object`) must now
  be named in full: `...` sits between them, and the constructors reject
  anything passed through it. In particular `data` must be named, as in
  `rc_nls(discharge, stage, data = thompson)`.

* Stage and discharge are named in full throughout, rather than `h` and `q`:
  * the constructors take `discharge` and `stage` (formerly `q` and `h`), as
    in `rc_nls(discharge, stage, data = thompson)`;
  * `predict()` and the `*_limits_2seg()` functions take `stage` (formerly
    `hpred`), and return it as the column `stage` (formerly `h`);
  * the gaugings stored on a fit are `gaugings`, with columns `discharge` and
    `stage` (formerly `qh_obs`, with `q` and `h`), and the fitted models'
    formulas use the same names;
  * `thompson` has columns `stage`, `discharge` and `uncertainty_pct`
    (formerly `h`, `q` and `uq`).

* `thompson$uncertainty_pct` is now documented as what it is: a
  percentage of the discharge at two standard deviations, as reported by the
  Water Survey of Canada. It was previously described as a discharge
  uncertainty, and the tests used `1 / uq^2` as weights as though it were a
  standard deviation in cubic meters per second.

* The `"sim"` interval method is removed, along with `sim_limits_2seg()`.
  Drawing parameters from their asymptotic normal distribution and pushing
  each draw through the model too often gave impossible curves (negative or
  astronomically large discharges) where a segment is poorly identified,
  making the limits meaningless. Use `method = "boot"` instead.

* By default, `rc_2seg_nls()` now tries 10 starting breakpoints spread across
  the search range and keeps the most likely fit, rather than starting once
  from the middle of the range. `kstart` may also be a vector of starting
  values to try. The fit is sensitive to where it starts: from some starts
  `nls()` fails, from others it stops at a local optimum. So fits that used
  to fail may now succeed, and fits that stopped at a local optimum may now
  find a better one, at a different breakpoint. What each start led to is
  recorded in `kstart_search`, and `settings$kstart` records `kstart` as
  supplied, so that `boot_limits_2seg()` repeats the same search for every
  resample. That costs one fit per start per resample, but shortcuts that
  start each resample only from the original estimate turned out to miss the
  resample's best fit too often, which would understate the uncertainty in
  the breakpoint.

* `rc_2seg_nls()` no longer takes `conflev` or `predlev`. They were stored on
  the fit but never used; give the levels to `predict()`.

* Fit objects are restructured, the same way for every model (see
  `?rating_curve`):
  * `pars` is a named list of the estimated curve parameters, one element per
    parameter type. For two-segment fits each element holds one value per
    segment, so `pars$b` is `c(b1, b2)`. `rc_poly()` names its coefficients
    `b0`, `b1`, ...; `rc_loess()` has none.
  * `settings` holds the arguments the fit was made with. It replaces
    `fit_args` and the separate `wts_code`, `formula` and `configuration`
    elements, and the loess tuning values formerly in `pars`.
  * `qh_obs` has columns `q` and `h` for every model (formerly `qobs` and
    `hobs` for the one-segment models).
  * The bias-corrected coefficients of the log-scale fits move from `pars` to
    `a_corrected`, and the `rc_gnls()` variance parameter from `pars$t_gnls`
    to `var_pars$power`.

* Package dependencies: R >= 4.0.0 is now required (was 4.1). tidyr and nls2
  are no longer used, and MASS is no longer suggested.

## New features

* `wts` is looked up in `data` like `discharge` and `stage`, so weights can
  be written in terms of its columns: `wts = 1 / uncertainty_sd^2`.

* `coef()` returns the estimated parameters as a flat named vector.

* Proportional weights (`wts_code = "prop"`) are fitted more robustly:
  * each reweighting round starts from the previous round's estimates, rather
    than from the initial starting values;
  * the rounds stop when no fitted discharge changes by more than `wts_tol`
    (relative), rather than when no coefficient does, which was unstable for
    an offset near zero;
  * running out of rounds (`wts_maxiter`) gives a warning, and every fit
    records the number of rounds and whether they converged in `irls`.
    `boot_limits_2seg()` treats a resample that does not converge as failed,
    and redraws it.

  Fits change only at the level of the tolerance: by at most about 2e-6,
  relative, on the Thompson and Sauze gaugings.

## Bug fixes

* `rc_poly(degree = 1)` fitted spurious quadratic and cubic terms.
* `rc_2seg_nls()` ignored `nls_maxiter`.
* User-supplied weights (`wts_code = "spec"`) fell out of step with the
  gaugings when any gauging had a missing stage or discharge.
* `rc_loess()` stored `NULL` for its residual standard error and degrees of
  freedom; it now stores `rse` and `enp`, the equivalent number of
  parameters.
* `print()` returns the fit invisibly, as documented.

# CSHShydRometry 0.0.1

* First packaged version. Derived from rating-curve functions written by Dan
  Moore, restructured into a model–predict form: every model is fitted by an
  `rc_*()` constructor and summarised by a `predict()` method returning the
  same columns, so results from different approaches stack with `rbind()`.
  Adds bootstrap and simulation interval methods for two-segment curves. The
  Bayesian methods, a breakpoint-averaged interval method and plotting
  helpers were left out.

* Scripts written against the model–predict functions can install this
  version with `remotes::install_github("CSHS-CWRA/CSHShydRometry@v0.0.1")`.
