# Design notes

Decisions about how CSHShydRometry is built, why, and what they leave open.
They are here so that a contributor can extend the package without
rediscovering the reasoning, and can tell a deliberate choice from an
accident. When a decision changes, update its entry rather than adding a new
one, and say in NEWS.md what changed for users.

## Weighting schemes

**Decision.** How the scatter of the gaugings is modelled is chosen by a
single `wts` argument taking a weighting object: `wts_none()`, `wts_prop()`,
`wts_power()` or `wts_spec()`.

- `wts_prop()` is scatter proportional to the fitted flow. It is a *known*
  variance function: nothing about the scatter is estimated beyond its
  scale. Every fitting function supports it, by reweighting in rounds (fit,
  recompute the weights from the fitted values, refit), in
  `reweight_in_rounds()`.
- `wts_power()` is scatter proportional to a power of the fitted flow, with
  the power *estimated*. Estimating a variance parameter needs a likelihood,
  which `nlme::gnls()` provides, so it is available in `rc_power()` only.
- Both are built on an internal `wts_nlme()`, which records the variance
  function in nlme's terms: `varPower(fixed = 1)` and `varPower()`. This is
  an internal convenience while `gnls()` does the estimating, not a
  commitment to nlme: see below.

**Why.**

- One argument, rather than a code plus settings that apply only to some
  codes, keeps the argument list short and avoids arguments whose meaning
  depends on another argument.
- A fixed variance function and an estimated one are different models, so
  they are different functions. An exponent argument that is sometimes fixed
  and sometimes estimated was tried and dropped.
- `wts_prop()` is fitted by our own rounds, not by `nlme::gnls()`, for two
  reasons. It has to work in every fitting function, and most are not fitted
  with nlme (`loess()`, `lm()`, and `nls()` with bounds on the parameters).
  And on small datasets `gnls()` was unreliable: on 14 gaugings it stopped in
  different places from different starts, and failed outright with a tight
  tolerance, where the rounds reached the same estimate every time. The nlme
  object *describes* the model; it does not fit it.

**What it leaves open.**

- Further estimated shapes for the scatter, such as an exponential
  (`wts_exp()`), would sit beside `wts_power()`. A public interface to
  nlme's variance functions is deliberately not planned, since the package
  may move away from nlme.
- `gnls()` is an implementation choice, not a statistical necessity. The
  power could be estimated for any fitting function, for example by
  maximising the profile likelihood over the power, refitting at each value
  with the rounds (`profile_loglik()` already exists for the breakpoint
  search). That would also extend `wts_power()` beyond `rc_power()`, and
  could replace `gnls()` altogether.
- Combining schemes, such as known gauging uncertainty plus scatter that
  grows with the flow (nlme's `varComb()`), would be a new function, such as
  `wts_comb()`.

## Naming

**Decision.**

- `rc_*()` functions fit models, and only they use the prefix. Internal
  helpers do not.
- Fitting functions are named for the model, not the fitting algorithm:
  `rc_power()`, not `rc_nls()`. A `_log` suffix marks a fit on the log-log
  scale.
- Multi-segment curves put the number of segments first: `rc_2seg_power()`,
  leaving room for `rc_2seg_*()` and `rc_pseg_*()`.
- Arguments and columns are named in full (`discharge`, `stage`), not with
  extreme abbreviations (`q`, `h`).

## Arguments

**Decision.**

- In exported functions, `...` sits straight after the required arguments,
  so every optional argument must be named in full. Fitting functions reject
  anything passed through `...` with `rlang::check_dots_empty()`.
- Arguments with one value per gauging (`discharge`, `stage`, and the values
  of `wts_spec()`) are evaluated in `data`, so they can refer to its columns.

## What a fit contains

**Decision.** Every fit has `gaugings`, `pars`, `settings`, `rse` and
`model`; see `?rating_curve`.

- `pars` holds the parameters of the *curve* only, as a named list with one
  element per parameter type and, for multi-segment curves, one value per
  segment. Elements may differ in length, where a configuration fixes a
  parameter. `coef()` flattens it.
- Parameters of the *scatter*, such as the power estimated under
  `wts_power()`, are not curve parameters: they live in the weighting
  scheme, `fit$wts`.
- `settings` records the arguments as given, so that a fit can be refitted
  (the bootstrap does this).

## Output of `predict()`

**Decision.** `predict()` returns a tibble with the same columns whatever the
model, method or weighting: `stage`, `fit`, and, when asked for,
`ci_lwr`/`ci_upr` and `pi_lwr`/`pi_upr`. A quantity that cannot be computed
is `NA`, not a missing column.

**Why.** Results from different approaches then stack with `rbind()`. Tables
are always tibbles, and the package imports from tibble so that loading it
loads tibble: otherwise a tibble behaves like a data frame (printing, `[`,
partial matching with `$`) until something else loads tibble, and the same
script can behave differently from one session to the next.

## Fits on the log-log scale

**Decision.** `rc_power_log()` takes no `wts`, and holds the stage of zero
flow fixed when `zero_flow_stage` is given.

**Why.** Equal scatter on the log scale already means scatter proportional to
the flow, which is usually the reason to weight. Taking logs does not always
even out the scatter, though, so weighting on the log scale is a possible
extension. A known `zero_flow_stage` makes the model linear in its other
parameters, fitted exactly by `lm()`. Fixing it at an *estimate* gives the
same curve with limits that leave out its uncertainty; that is how the former
`rc_log_ols()` behaved, and it is reproduced by fitting twice rather than
built in.

## Two-segment fits

**Decision.**

- `rc_2seg_power()` tries several starting breakpoints and keeps the most
  likely fit, because the fit is sensitive to where it starts.
- `predict()` offers `method = "delta"` (the default, fast, unreliable near
  the breakpoint) and `method = "boot"`.
- The bootstrap repeats the full breakpoint search for every resample.

**Why.**

- Starting each resample from the original breakpoint, or from the optima
  the original search found, was tried. It missed the resample's best fit in
  up to half of the resamples, which would understate the uncertainty in the
  breakpoint.
- A method drawing parameters from their asymptotic normal distribution
  (`"sim"`) was tried and removed: where a segment is poorly identified, the
  draws too often give impossible curves.
