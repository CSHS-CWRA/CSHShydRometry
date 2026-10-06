# Fit rating curve using loess smoother

Fit rating curve using loess smoother

## Usage

``` r
rc_loess(
  discharge,
  stage,
  ...,
  data = NULL,
  degree = 2,
  span = 0.75,
  extrapolate = TRUE,
  variance = var_none()
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

- degree:

  Degree of local polynomials (1 or 2).

- span:

  Smoothing parameter, passed to
  [`stats::loess()`](https://rdrr.io/r/stats/loess.html).

- extrapolate:

  Single logical; allow extrapolation beyond the observed stage range?
  Default is `TRUE`. If `FALSE`, predictions outside the range of the
  gaugings are `NA`. (This sets loess's `surface` to `"direct"` or
  `"interpolate"`; see
  [`stats::predict.loess()`](https://rdrr.io/r/stats/predict.loess.html).)

- variance:

  **\[experimental\]** How the variance of the gaugings about the curve
  is modelled:
  [`var_none()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  (or `"none"`, the default),
  [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  (or `"prop"`), or
  [`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  with the variances. See
  [variance](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md).
  Under
  [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  the loess curve is refitted in rounds, like the parametric fits.

## Value

An `rc_loess` object; see
[rating_curve](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md)
for its contents. A loess curve has no parameters, so `curve_parameters`
is an empty list.

## Examples

``` r
fit <- rc_loess(discharge, stage, data = thompson)
predict(fit, new_stage = c(1, 3, 6), conflev = 0.95)
#> # A tibble: 3 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  250.   228.   273.
#> 2     3  823.   796.   850.
#> 3     6 2216.  2192.  2240.
```
