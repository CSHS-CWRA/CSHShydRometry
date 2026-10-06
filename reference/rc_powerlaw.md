# Fit a power-law rating curve on the original scale

Fits the power law \\Q = a (h - c)^b\\ by nonlinear least squares on the
original scale of discharge, with the scatter of the gaugings modelled
as set by `variance`. The curve estimates the mean discharge at each
stage.

## Usage

``` r
rc_powerlaw(
  discharge,
  stage,
  ...,
  data = NULL,
  variance = var_none(),
  offset = NULL,
  control = stats::nls.control(maxiter = 1000, tol = 1e-06)
)
```

## Arguments

- discharge:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Discharge: a vector, or an expression evaluated in `data`, such as a
  column name.

- stage:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  Stage: a vector, or an expression evaluated in `data`, such as a
  column name.

- ...:

  Must be empty. Present so that every argument after it has to be named
  in full.

- data:

  Optional data frame in which `discharge` and `stage` are evaluated.

- variance:

  **\[experimental\]** How the variance of the gaugings about the curve
  is modelled:
  [`var_none()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  (or `"none"`, the default),
  [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  (or `"prop"`),
  [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  (or `"power"`), or
  [`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  with the variances. See
  [variance](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md).

- offset:

  The offset, \\c\\ in the formula: `NULL` (the default) to estimate it,
  or a known value, below every gauged stage, to hold it fixed.

- control:

  Settings for [`stats::nls()`](https://rdrr.io/r/stats/nls.html), as
  from
  [`stats::nls.control()`](https://rdrr.io/r/stats/nls.control.html).

## Value

An `rc_powerlaw` object; see
[rating_curve](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md)
for its contents.

## The offset

\\c\\ is the offset: for a single power law, the stage at which the flow
would stop. It is estimated unless it is given, through `offset`, say
from a survey of the control. A given offset is held fixed, so the
limits from [`predict()`](https://rdrr.io/r/stats/predict.html) carry no
uncertainty in it. Fixing it at the value estimated by a first fit
treats an estimate as known: the curve is the same, but the limits are
narrower than they should be.

## Examples

``` r
fit <- rc_powerlaw(discharge, stage, data = thompson)
predict(fit, new_stage = c(1, 3, 6), conflev = 0.95)
#> # A tibble: 3 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  243.   222.   264.
#> 2     3  832.   814.   850.
#> 3     6 2209.  2188.  2230.

# constant coefficient of variation instead of constant variance
fit_prop <- rc_powerlaw(discharge, stage, data = thompson, variance = "prop")
predict(fit_prop, new_stage = c(1, 3, 6), conflev = 0.95)
#> # A tibble: 3 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  252.   249.   255.
#> 2     3  816.   806.   826.
#> 3     6 2205.  2178.  2232.

# the offset known, say from a survey of the control
rc_powerlaw(discharge, stage, data = thompson, variance = "prop", offset = -1.3)
#> Rating curve model.
#> - Method: rc_powerlaw 
```
