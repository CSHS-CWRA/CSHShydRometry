# Estimated parameters of a rating curve

The contents of `object$curve_parameters` as a flat named numeric
vector. For a single-segment curve the names are those of
`curve_parameters`; for a multi-segment curve they carry the segment
number, as in the model formula (`a1`, `b2`, `k`, ...).

## Usage

``` r
# S3 method for class 'rating_curve'
coef(object, ...)

# S3 method for class 'rc_2seg_powerlaw'
coef(object, ...)
```

## Arguments

- object:

  A rating curve fit.

- ...:

  Ignored.

## Value

A named numeric vector; empty for
[`rc_loess()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_loess.md).

## Examples

``` r
coef(rc_powerlaw(discharge, stage, data = thompson))
#>          a          b          c 
#> 80.6686509  1.7131161 -0.9044142 
```
