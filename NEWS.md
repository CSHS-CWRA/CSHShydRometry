# CSHShydRometry (development version)

## Breaking changes

* Fitting functions are named for the model they fit, not the algorithm
  that fits it, with the number of segments first for multi-segment curves:
  * `rc_nls()` is now `rc_powerlaw()`: a power law, by least squares on the
    original scale;
  * `rc_log_nls()` is now `rc_powerlaw_log()`: a power law, by least squares on
    the log-log scale;
  * `rc_nls_2seg()` is now `rc_2seg_powerlaw()`.

  Their classes are renamed to match. This leaves room for more two-segment
  methods (`rc_2seg_*()`) and for curves with more segments.

* `rc_gnls()` is removed. The scatter it modelled, proportional to the flow
  raised to an estimated power, is now a variance scheme of its own,
  `var_power()`: `rc_powerlaw(variance = var_power())`, or `variance = "power"`. It adds
  a parameter, the power, so it is fitted by generalised least squares with
  `nlme::gnls()` as before, and is available in `rc_powerlaw()` only. The fit
  now starts from the fit under `var_prop()`, rather than from a straight
  line on the log-log scale, and records the estimated power in
  `fit$variance$exponent`. `rc_gnls()` also took any nlme variance function
  (`var_type`); only the power of the mean, which every example used, is
  kept.

* Limits for fits under `var_power()` are now computed by the delta method,
  as for the other power-law fits, rather than simulated by
  `nlraa::predict_gnls()`. They are the same every time `predict()` is
  called, and treat the estimated power as known. The simulated confidence
  limits were also too wide: the parameters of a power law are strongly
  correlated, and draws from their approximate normal distribution stray
  where the curve is far from linear in them. On the Thompson gaugings, a
  case bootstrap agrees with the delta method. nlraa is no longer a
  dependency, and `predict()` no longer passes `...` on to it.

* `rc_log_ols()` is removed. It gave the same estimates as `rc_log_nls()`,
  but treated the estimated `c` as known when computing limits, which left
  its uncertainty out. `rc_powerlaw_log()` instead gains an `offset`
  argument for `c`: `NULL` (the default) estimates it, and a value holds it
  fixed, for an offset known from a survey. Fixing it at the
  estimate of `c` reproduces the old behaviour. `rc_powerlaw()` takes
  `offset` too.

* `c` is called the offset, not the stage of zero flow. The two coincide for
  a single power law, but not for a segment of a multi-segment curve.

* In every `predict()` method, and in `delta_limits_2seg()` and
  `boot_limits_2seg()`, the stages to predict at are now given as
  `new_stage`, not `stage`.

* `predict()` methods no longer pass `...` on to the functions they call.
  An argument a method does not take, such as `method` for a
  single-segment fit or a misspelt `conflev`, is now an error rather than
  silently ignored.

* How the gaugings scatter about the curve is chosen by a single argument,
  `variance`, in `rc_powerlaw()`, `rc_poly()`, `rc_loess()` and
  `rc_2seg_powerlaw()`. It replaces `wts_code`, `wts`, `wts_tol` and
  `wts_maxiter`, whose meanings depended on one another. `variance` takes a
  variance scheme: `var_none()` (the default), `var_prop(tol, maxiter)`,
  `var_power()` or `var_spec(values)`, or the shorthand `"none"`, `"prop"`
  or `"power"`; see `?variance`. The schemes describe variances, not
  weights: the fits weight each gauging by the reciprocal of its variance.
  So `var_spec()` takes the variances themselves, such as
  `var_spec(sauze$uncertainty_sd^2)`, as an ordinary vector. Describing
  variances leaves room to combine sources of scatter, such as measurement
  error and the rest, by adding them. This interface is experimental.

* New `fitted()` and `residuals()` methods for every fit. `residuals()` gives
  observed minus fitted discharge, or, with `type = "pearson"`, residuals
  scaled by their modelled standard deviation, for checking the variance
  scheme and normality. Use these rather than the `model` element, which
  depends on how a fit is computed and may change.

* The settings passed to `nls()` are given as a single `control` argument,
  as from `stats::nls.control()`, like `nls()` itself. It replaces `nls_tol`
  and `nls_maxiter` in `rc_powerlaw()` and `rc_2seg_powerlaw()`, and `tol` in
  `rc_powerlaw_log()`.

* `rc_2seg_powerlaw()`'s `config` argument is renamed `combine`, since it says
  how the two segments combine above the breakpoint, and its values
  `"piecewise"` and `"compound"` are renamed `"replace"` (the upper power law
  takes over from the lower) and `"add"` (it adds to the lower). Both kinds
  of curve are piecewise, so the old values did not tell them apart.

* Every argument after the mandatory ones (`discharge` and `stage`, or `object`) must now
  be named in full: `...` sits between them, and the constructors reject
  anything passed through it. In particular `data` must be named, as in
  `rc_powerlaw(discharge, stage, data = thompson)`.

* Stage and discharge are named in full throughout, rather than `h` and `q`:
  * the constructors take `discharge` and `stage` (formerly `q` and `h`), as
    in `rc_powerlaw(discharge, stage, data = thompson)`;
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

* By default, `rc_2seg_powerlaw()` now tries 10 starting breakpoints spread across
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

* `rc_2seg_powerlaw()` no longer takes `contcons`, which chose the parameter
  carrying the continuity constraint. Only `"a"` was implemented, and that
  is what `combine = "replace"` does.

* `rc_2seg_powerlaw()` no longer takes `conflev` or `predlev`. They were stored on
  the fit but never used; give the levels to `predict()`.

* Fit objects are restructured, the same way for every model (see
  `?rating_curve`):
  * `pars` is a named list of the estimated curve parameters, one element per
    parameter type. For two-segment fits each element holds one value per
    segment, so `pars$b` is `c(b1, b2)`. `rc_poly()` names its coefficients
    `b0`, `b1`, ...; `rc_loess()` has none.
  * `weights_used` holds the weights the final model was fitted with
    (formerly `weights`).
  * `settings` holds the arguments the fit was made with. It replaces
    `fit_args` and the separate `wts_code`, `formula` and `configuration`
    elements, and the loess tuning values formerly in `pars`.
  * `gaugings` has the same columns, `discharge` and `stage`, for every model
    (formerly `qh_obs`, with `qobs` and `hobs` for the one-segment models and
    `q` and `h` for the two-segment one).
  * The bias-corrected coefficients of the log-scale fits move from `pars` to
    `a_corrected`, and the power formerly in `pars$t_gnls` of `rc_gnls()`
    fits to `variance$exponent` (see `rc_gnls()` above).

* Package dependencies: R >= 4.0.0 is now required (was 4.1). tidyr and nls2
  are no longer used, and MASS is no longer suggested.

## New features

* Two vignettes: `vignette("fitting")` walks through the fitting functions,
  from choosing gaugings and looking at their scatter to two-segment curves,
  and `vignette("uncertainty")` covers confidence and prediction limits, how
  the weighting shapes them, and where they mislead.

* `thompson` gains `rating_table`, the Water Survey's rating table in force
  at each gauging, where recorded (2020 onwards). The rating has shifted
  over the decades, and its documentation now says so.

* Tables are always tibbles: `predict()` and the limits functions return
  them, fits store their gaugings as one, and `thompson` is one. tibble is
  now a dependency (it was suggested). Previously the output was a tibble
  only when tibble was installed, so the same script could behave
  differently from one machine to another; for example, `$` partial-matches
  column names on a data frame but not on a tibble.

* `coef()` returns the estimated parameters as a flat named vector.

* Proportional weights (`var_prop()`) are fitted more robustly:
  * each reweighting round starts from the previous round's estimates, rather
    than from the initial starting values;
  * the rounds stop when no fitted discharge changes by more than `tol`
    (relative), rather than when no coefficient does, which was unstable for
    an offset near zero;
  * running out of rounds (`maxiter`) gives a warning, and every fit
    records the number of rounds and whether they converged in `irls`.
    `boot_limits_2seg()` treats a resample that does not converge as failed,
    and redraws it;
  * `rc_loess()` now reweights in rounds too. It used to reweight once, from
    an unweighted fit, so its weights did not match its own fitted values.

  Fits change only at the level of the tolerance: by at most about 2e-6,
  relative, on the Thompson and Sauze gaugings.

## Bug fixes

* `rc_poly(degree = 1)` fitted spurious quadratic and cubic terms.
* `rc_2seg_powerlaw()` ignored `nls_maxiter`.
* `nls_tol` had no effect in `rc_2seg_powerlaw()`: its `"port"` algorithm ignores
  `tol`. `control` now says so, and passes on the port algorithm's own
  settings, such as `rel.tol`.
* User-supplied weights (now `var_spec()`) fell out of step with the
  gaugings when any gauging had a missing stage or discharge.
* `rc_loess()` stored `NULL` for its residual standard error and degrees of
  freedom; it now stores `rse` and `enp`, the equivalent number of
  parameters.
* `print()` returns the fit invisibly, as documented.

# CSHShydRometry 0.0.1

* First packaged version. Derived from rating-curve functions written by R. Dan
  Moore, restructured into a model–predict form: every model is fitted by an
  `rc_*()` constructor and summarised by a `predict()` method returning the
  same columns, so results from different approaches stack with `rbind()`.
  Adds bootstrap and simulation interval methods for two-segment curves. The
  Bayesian methods, a breakpoint-averaged interval method and plotting
  helpers were left out.

* Scripts written against the model–predict functions can install this
  version with `remotes::install_github("CSHS-CWRA/CSHShydRometry@v0.0.1")`.
