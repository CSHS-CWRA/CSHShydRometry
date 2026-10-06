# Variance schemes for rating-curve fits

How the variance of the gaugings about the curve is modelled. Pass one
of these as the `variance` argument of
[`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md),
[`rc_poly()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_poly.md),
[`rc_loess()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_loess.md)
or
[`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md).
The fits weight each gauging by the reciprocal of its variance.

## Usage

``` r
var_none()

var_prop(tol = 1e-06, maxiter = 100)

var_power()

var_spec(values)

# S3 method for class 'rc_variance'
print(x, ...)
```

## Arguments

- tol:

  Convergence tolerance: the reweighting stops once no fitted discharge
  changes by more than this fraction from one round to the next.

- maxiter:

  Maximum number of reweighting rounds.

- values:

  The variances, one per gauging, as a numeric vector: for example
  `sd^2`, where `sd` is the standard uncertainty of each discharge. They
  are evaluated where `var_spec()` is called, as an ordinary argument,
  not looked up in the fit's `data`; write `sauze$uncertainty_sd^2`, not
  `uncertainty_sd^2`. The variances of gaugings that the fit drops, for
  a missing stage or discharge, are dropped with them.

- x:

  An `"rc_variance"` object.

- ...:

  Ignored.

## Value

An object of class `"rc_variance"`.

## Details

**\[experimental\]** The interface for modelling the scatter is likely
to grow, and may change.

- `var_none()`: the variance is the same at every flow (ordinary least
  squares).

- `var_prop()`: the standard deviation is proportional to the flow, that
  is, a constant coefficient of variation. The variances, proportional
  to `fitted^2`, depend on the fit itself, so the fit is repeated in
  rounds: fit, recompute the weights from the fitted values, refit, each
  round starting from the previous round's estimates. The rounds stop
  once no fitted discharge changes by more than a fraction `tol` from
  one round to the next, or after `maxiter` rounds, with a warning.

- `var_power()`: the standard deviation is proportional to a power of
  the flow, with the power estimated along with the curve. This adds a
  parameter, so the fit is made by generalised least squares with
  [`nlme::gnls()`](https://rdrr.io/pkg/nlme/man/gnls.html), starting
  from the fit under `var_prop()`. It is available in
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md)
  only: that is a matter of implementation, as the other fitting
  functions are not made with
  [`nlme::gnls()`](https://rdrr.io/pkg/nlme/man/gnls.html). The estimate
  can be unstable with few gaugings; compare it with the fit under
  `var_prop()`. The estimated power is kept on the fit, in
  `fit$variance$exponent`. The limits from
  [`predict()`](https://rdrr.io/r/stats/predict.html) treat it as known,
  so they leave out the uncertainty in the power.

- `var_spec()`: the relative variances of the gaugings are known,
  typically from their reported uncertainties, and given directly. Only
  their relative sizes matter: the fit estimates a common scale factor,
  the residual standard error, which is 1 when the gaugings scatter
  exactly as their variances say. A new gauging has no variance given,
  so prediction limits are returned as `NA`.

The strings `"none"`, `"prop"` and `"power"` are shorthand for
`var_none()`, `var_prop()` and `var_power()` with their defaults.

## Examples

``` r
rc_powerlaw(discharge, stage, data = thompson, variance = var_prop())
#> Rating curve model.
#> - Method: rc_powerlaw 

# the same, with the defaults
rc_powerlaw(discharge, stage, data = thompson, variance = "prop")
#> Rating curve model.
#> - Method: rc_powerlaw 

# the power of the flow estimated too
fit <- rc_powerlaw(discharge, stage, data = thompson, variance = var_power())
fit$variance$exponent
#> [1] 0.999845

# variances from each gauging's reported uncertainty, leaving out the one
# reported as 0.0226%, probably a fraction entered as a percentage (see
# ?thompson): it would carry over 99.9% of the total weight
d <- thompson[!is.na(thompson$uncertainty_pct) & thompson$uncertainty_pct > 1, ]
d$uncertainty_sd <- d$uncertainty_pct / 100 * d$discharge / 2
rc_powerlaw(discharge, stage, data = d, variance = var_spec(d$uncertainty_sd^2))
#> Rating curve model.
#> - Method: rc_powerlaw 
```
