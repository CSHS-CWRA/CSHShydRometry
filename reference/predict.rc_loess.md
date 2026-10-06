# Predict method for rc_loess objects

Predict method for rc_loess objects

## Usage

``` r
# S3 method for class 'rc_loess'
predict(object, ..., new_stage = NULL, conflev = NULL, predlev = NULL)
```

## Arguments

- object:

  An rc_loess object.

- ...:

  Must be empty. Present so that every argument after it has to be named
  in full; a misspelt or unsupported argument is an error, not silently
  ignored.

- new_stage:

  Stages at which to predict discharge. If `NULL`, a grid of 1000
  equally spaced stages spanning the observed range is used.

- conflev:

  The confidence level for the confidence limits; if NULL, no confidence
  limits are returned.

- predlev:

  The prediction level for the prediction limits; if NULL, no prediction
  limits are returned.

## What the limits assume

The limits use the t distribution, so they assume the scatter of the
gaugings about the curve is normal (on the log scale, for
[`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md)),
with the spread the weighting scheme describes. Confidence limits depend
on this only mildly, because estimates average over many gaugings.
Prediction limits depend on it directly: if the scatter is skewed or has
heavy tails, they can miss, especially at high levels such as 0.99.
Check the residuals before relying on them.

Even when the scatter is normal, most limits are approximations: a 95%
interval covers the truth roughly, not exactly, 95% of the time. That is
because the curve is nonlinear in its parameters and is approximated by
a straight line about the estimates (the delta method), because the
weights are themselves estimated from the fit (under
[`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
and
[`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)),
or, for
[`rc_loess()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_loess.md),
because of the smoother's approximations. The approximation is good with
plenty of gaugings, and poorer with few, or beyond the range of the
gaugings. The limits are exact only for fits that are linear in their
parameters with weights fixed in advance:
[`rc_poly()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_poly.md)
with
[`var_none()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
or
[`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
and
[`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md)
with `offset` given.
