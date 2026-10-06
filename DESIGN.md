# Design notes

Decisions about how CSHShydRometry is built, why, and what they leave
open. They are here so that a contributor can extend the package without
rediscovering the reasoning, and can tell a deliberate choice from an
accident. When a decision changes, update its entry rather than adding a
new one, and say in NEWS.md what changed for users.

## Variance schemes

**Decision.** How the scatter of the gaugings is modelled is chosen by a
single `variance` argument taking a variance scheme:
[`var_none()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
[`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
[`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
or
[`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md).
The schemes describe *variances*; the fits weight each gauging by the
reciprocal of its variance, and users never handle weights. The
interface is marked experimental.

- [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  is a standard deviation proportional to the fitted flow. It is a
  *known* variance function: nothing about the scatter is estimated
  beyond its scale. Every fitting function supports it, by reweighting
  in rounds (fit, recompute the weights from the fitted values, refit),
  in `reweight_in_rounds()`.
- [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  is a standard deviation proportional to a power of the fitted flow,
  with the power *estimated*. Estimating a variance parameter needs a
  likelihood, which
  [`nlme::gnls()`](https://rdrr.io/pkg/nlme/man/gnls.html) provides, so
  it is available in
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md)
  only.
- Both are built on an internal `var_nlme()`, which records the variance
  function in nlme’s terms: `varPower(fixed = 1)` and `varPower()`. This
  is an internal convenience while `gnls()` does the estimating, not a
  commitment to nlme: see below.

**Why.**

- Variance, not weight, is the quantity people think in: how much a
  gauging scatters. Weights were an artefact of least squares, and would
  not carry over to other ways of fitting.
- Variances add, and weights do not. Describing variances leaves room to
  combine sources of scatter later, with `+`, for example measurement
  error and the rest. An earlier interface, `wts`, took weights; it was
  renamed before the first release for this reason.
- One argument, rather than a code plus settings that apply only to some
  codes, keeps the argument list short and avoids arguments whose
  meaning depends on another argument.
- A fixed variance function and an estimated one are different models,
  so they are different functions. An exponent argument that is
  sometimes fixed and sometimes estimated was tried and dropped.
- [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  is fitted by our own rounds, not by
  [`nlme::gnls()`](https://rdrr.io/pkg/nlme/man/gnls.html), for two
  reasons. It has to work in every fitting function, and most are not
  fitted with nlme ([`loess()`](https://rdrr.io/r/stats/loess.html),
  [`lm()`](https://rdrr.io/r/stats/lm.html), and
  [`nls()`](https://rdrr.io/r/stats/nls.html) with bounds on the
  parameters). And on small datasets `gnls()` was unreliable: on 14
  gaugings it stopped in different places from different starts, and
  failed outright with a tight tolerance, where the rounds reached the
  same estimate every time. The nlme object *describes* the model; it
  does not fit it.

**What it leaves open.**

- Further estimated shapes for the scatter, such as an exponential
  (`var_exp()`), would sit beside
  [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md).
  A public interface to nlme’s variance functions is deliberately not
  planned, since the package may move away from nlme.
- `gnls()` is an implementation choice, not a statistical necessity. The
  power could be estimated for any fitting function, for example by
  maximising the profile likelihood over the power, refitting at each
  value with the rounds. That would also extend
  [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  beyond
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md),
  and could replace `gnls()` altogether.
- Combining schemes, such as known gauging uncertainty plus scatter that
  grows with the flow, by adding them with `+`. Once two components each
  have a scale to estimate, weighted least squares cannot fit them, and
  the fit moves to likelihood, with least squares as the special case.

## Naming

**Decision.**

- `rc_*()` functions fit models, and only they use the prefix. Internal
  helpers do not.
- Fitting functions are named for the model, not the fitting algorithm:
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md),
  not `rc_nls()`. A `_log` suffix marks a fit on the log-log scale.
- Multi-segment curves put the number of segments first:
  [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md),
  leaving room for `rc_2seg_*()` and `rc_pseg_*()`.
- Arguments and columns are named in full (`discharge`, `stage`), not
  with extreme abbreviations (`q`, `h`).
- In formulas, a power law is `Q = a (h - c)^b`, and the parameters are
  named `a`, `b` and `c` in `curve_parameters` and
  [`coef()`](https://rdrr.io/r/stats/coef.html). As an argument, though,
  `c` is named `offset`.

**Why.**

- `c` lines up with `a` and `b`. The alternative, `h_0`, reads as the
  stage of zero flow, which `c` is only for a single power law: for a
  segment of a multi-segment curve it is not.
- An argument `c =` reads as R’s [`c()`](https://rdrr.io/r/base/c.html),
  so the argument has a name in words. “Offset” holds for any segment;
  “zero flow stage” does not.

## Arguments

**Decision.**

- In exported functions, `...` sits straight after the required
  arguments, so every optional argument must be named in full. Fitting
  functions reject anything passed through `...` with
  [`rlang::check_dots_empty()`](https://rlang.r-lib.org/reference/check_dots_empty.html).
- `discharge` and `stage` are evaluated in `data`, so they can refer to
  its columns.
- In [`predict()`](https://rdrr.io/r/stats/predict.html), the stages to
  predict at are `new_stage`, distinct from the `stage` the curve was
  fitted to.
- [`predict()`](https://rdrr.io/r/stats/predict.html) methods check that
  `...` is empty too, rather than passing it on. A single-segment fit
  has no `method`, so `predict(fit, method = "boot")` is an error rather
  than delta-method limits that look like bootstrap ones. `method`
  arrives with a single-segment bootstrap. The two-segment
  [`predict()`](https://rdrr.io/r/stats/predict.html) passes `...` on to
  its limits function, for `B` and `seed`, and that function checks it.
- The variances given to
  [`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  are *not*: they are an ordinary vector, evaluated where
  [`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  is called. A variance scheme is an object in its own right, which can
  be made in one place and used in another, so an expression captured
  inside it would take its meaning from whichever fit’s `data` it later
  met: a column could silently stand in for a variable of the same name,
  and a misspelt column would only fail inside the fit. The convenience
  (not writing `data$`) is small. Adding evaluation in `data` later
  would not break anyone’s code, whereas removing it would, so it was
  left out until it is asked for.

## What a fit contains

**Decision.** Every fit has `gaugings`, `curve_parameters` and
`settings`, and, where it takes a `variance` argument, `variance`; see
[`?rating_curve`](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md).
Other elements, such as `model`, `rse`, `weights_used` and `irls`,
depend on how the fit is computed, and are documented as liable to
change. Users get fitted values, residuals and parameters through
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`coef()`](https://rdrr.io/r/stats/coef.html).

- `curve_parameters` holds the parameters of the *curve* only, as a
  named list with one element per parameter type and, for multi-segment
  curves, one value per segment. Elements may differ in length, where a
  configuration fixes a parameter.
  [`coef()`](https://rdrr.io/r/stats/coef.html) flattens it.
- Parameters of the *scatter*, such as the power estimated under
  [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
  are not curve parameters: they live in the variance scheme,
  `fit$variance`.
- `model`, the underlying model object, is kept out of the interface
  because the engine is expected to change, for example to likelihood
  fitting when variance schemes are combined.
- `settings` records the arguments as given, so that a fit can be
  refitted (the bootstrap does this).

## Output of `predict()`

**Decision.** [`predict()`](https://rdrr.io/r/stats/predict.html)
returns a tibble with the same columns whatever the model, method or
variance scheme: `stage`, `fit`, and, when asked for, `ci_lwr`/`ci_upr`
and `pi_lwr`/`pi_upr`. A quantity that cannot be computed is `NA`, not a
missing column.

**Why.** Results from different approaches then stack with
[`rbind()`](https://rdrr.io/r/base/cbind.html). Tables are always
tibbles, and the package imports from tibble so that loading it loads
tibble: otherwise a tibble behaves like a data frame (printing, `[`,
partial matching with `$`) until something else loads tibble, and the
same script can behave differently from one session to the next.

## Fits on the log-log scale

**Decision.**
[`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md)
takes no `variance`, and holds the offset fixed when `offset` is given.

**Why.** Equal scatter on the log scale already means scatter
proportional to the flow, which is usually the reason to model the
variance. Taking logs does not always even out the scatter, though, so a
variance scheme on the log scale is a possible extension. A known
`offset` makes the model linear in its other parameters, fitted exactly
by [`lm()`](https://rdrr.io/r/stats/lm.html). Fixing it at an *estimate*
gives the same curve with limits that leave out its uncertainty; that is
how the former `rc_log_ols()` behaved, and it is reproduced by fitting
twice rather than built in.

## Two-segment fits

**Decision.**

- [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md)
  tries several starting breakpoints and keeps the most likely fit,
  because the fit is sensitive to where it starts.
- [`predict()`](https://rdrr.io/r/stats/predict.html) offers
  `method = "delta"` (the default, fast, unreliable near the breakpoint)
  and `method = "boot"`.
- The bootstrap repeats the full breakpoint search for every resample.
- Of the fits from different starts, the one kept has the smallest loss:
  the quantity the fit minimises. Under
  [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
  that is the Gamma quasi-likelihood (negated), not the sum of squared
  relative residuals: reweighting in rounds settles where the
  quasi-likelihood is maximised.
- How the segments combine above the breakpoint is `combine = "replace"`
  or `"add"`.

**Why.**

- Starting each resample from the original breakpoint, or from the
  optima the original search found, was tried. It missed the resample’s
  best fit in up to half of the resamples, which would understate the
  uncertainty in the breakpoint.
- A method drawing parameters from their asymptotic normal distribution
  (`"sim"`) was tried and removed: where a segment is poorly identified,
  the draws too often give impossible curves.
- `combine` was earlier `config`, after BaRatin’s configuration matrix,
  which says which segments carry flow in which range of stage: a row
  per range, a column per segment. With two segments only two matrices
  make sense, `rbind(c(1, 0), c(0, 1))` (“replace”) and
  `rbind(c(1, 0), c(1, 1))` (“add”), so they are named in words.

**What it leaves open.**

- With more segments, two words will not be enough. `combine` could also
  accept a configuration matrix, with `"replace"` and `"add"` kept as
  shorthands, so no existing code breaks. A matrix is hard to write and
  read, though, so how users specify one needs more thought.
