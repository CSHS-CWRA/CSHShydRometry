# CSHShydRometry: Statistical Methods for Rating Curves

Fits stage-discharge rating curves by a range of statistical methods,
and attaches confidence and prediction limits to the fitted curve.
Single-segment curves may be fitted as a power law (by ordinary or
nonlinear least squares, on the log-log or the natural scale, or by
generalised nonlinear least squares), as a polynomial, or by loess.
Two-segment power-law curves are also supported: above an estimated
breakpoint, the upper segment either takes over from the lower one or
adds to it. Every fit is given a 'predict()' method returning the same
columns, so results from different approaches can be combined directly.

## Fitting a rating curve

Every model is fitted by a `rc_*()` constructor and evaluated by a
[`predict()`](https://rdrr.io/r/stats/predict.html) method. The
constructors take discharge and stage — as vectors, or as column names
with a `data` argument — and return an object carrying `"rating_curve"`
as its top-level class; see
[rating_curve](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md).

Single segment:

- [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md)
  — power law fitted by least squares on the original scale.

- [`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md)
  — power law fitted by least squares on the log-log scale.

In both, the offset \\c\\ is estimated, or held at a value given as
`offset`.

- [`rc_poly()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_poly.md),
  [`rc_loess()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_loess.md)
  — polynomial and loess alternatives.

Two segments, joined at an estimated breakpoint:

- [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md)
  — with `combine = "replace"`: the upper power law takes over from the
  lower at the breakpoint; with `combine = "add"`: it adds to the
  discharge carried at the breakpoint.

## Variance

The constructors that take `variance` offer four models of how the
gaugings scatter about the curve. **\[experimental\]** This interface
may change:

- [`var_none()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  — the same scatter at every flow (constant variance).

- [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  — scatter proportional to the flow.

- [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  — scatter proportional to a power of the flow, with the power
  estimated;
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md)
  only at this time.

- [`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  — variances supplied by the user, typically from reported gauging
  uncertainties. A new observation's scatter is then not identified by
  the fit, so prediction limits are returned as `NA`.

See
[variance](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md).

## Confidence and prediction limits

[`predict()`](https://rdrr.io/r/stats/predict.html) returns the same
columns whatever the model, the method or the weighting: `stage`, `fit`,
and — when `conflev` or `predlev` is given — `ci_lwr`/`ci_upr` and
`pi_lwr`/`pi_upr`. Quantities that cannot be computed come back as `NA`
rather than as missing columns, so results from different approaches
stack directly.

For two-segment curves,
[`predict.rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/predict.rc_2seg_powerlaw.md)
takes a `method`. The default, `"delta"`, is fast but unreliable near
the breakpoint, where the mean function is not differentiable; `"boot"`
costs a refit per resample and behaves much better there. See
[`predict.rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/predict.rc_2seg_powerlaw.md)
for the detail.

## Data

[thompson](https://cshs-cwra.github.io/CSHShydRometry/reference/thompson.md)
ships with the package and suits the single-segment models. The
two-segment examples prefer the Ardeche at Sauze,
[`RBaM::SauzeGaugings`](https://rdrr.io/pkg/RBaM/man/SauzeGaugings.html),
which has a clearer change of control; RBaM is a suggested dependency
used only as a source of that data.

## See also

Useful links:

- <https://github.com/CSHS-CWRA/CSHShydRometry>

- Report bugs at <https://github.com/CSHS-CWRA/CSHShydRometry/issues>

## Author

Vincenzo Coia (maintainer, <vincenzo.coia@gmail.com>), R. Dan Moore, and
Paul Whitfield.
