# Fit polynomial rating curve

Fit polynomial rating curve

## Usage

``` r
rc_poly(discharge, stage, ..., data = NULL, degree = 2, variance = var_none())
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

  Polynomial degree; positive whole number. Default 2.

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

## Value

An `rc_poly` object; see
[rating_curve](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md)
for its contents. The coefficients are named `b0`, `b1`, ... by power of
`stage`.

## Examples

``` r
fit <- rc_poly(discharge, stage, data = thompson, degree = 2)
predict(fit, new_stage = c(1, 3, 6), conflev = 0.95)
#> # A tibble: 3 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  244.   222.   266.
#> 2     3  835.   816.   855.
#> 3     6 2201.  2181.  2221.
```
