# Fit a power-law rating curve on the log-log scale

Fits the power law \\Q = a (h - c)^b\\ by least squares on the log-log
scale, `log(discharge) ~ log(a) + b * log(stage - c)`. Equal scatter on
that scale means the scatter in discharge is proportional to the flow,
and the back-transformed curve estimates the geometric mean discharge at
each stage (also the median, when the scatter on the log scale is
symmetric).

## Usage

``` r
rc_powerlaw_log(
  discharge,
  stage,
  ...,
  data = NULL,
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

- offset:

  The offset, \\c\\ in the formula: `NULL` (the default) to estimate it,
  or a known value, below every gauged stage, to hold it fixed.

- control:

  Settings for [`stats::nls()`](https://rdrr.io/r/stats/nls.html), as
  from
  [`stats::nls.control()`](https://rdrr.io/r/stats/nls.control.html).
  Used only when the offset is estimated.

## Value

An `rc_powerlaw_log` object; see
[rating_curve](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md)
for its contents. The back-transformed coefficient `a` estimates the
geometric mean discharge; `a_corrected` holds two bias-corrected
versions for the arithmetic mean: `nbc`, assuming lognormal errors, and
`dbc`, Duan's smearing estimate.

## The offset

\\c\\ is the offset: for a single power law, the stage at which the flow
would stop. When it is estimated, the fit is nonlinear in \\c\\ and is
made with [`stats::nls()`](https://rdrr.io/r/stats/nls.html). When it is
given, through `offset`, the model is linear in its remaining parameters
on the log-log scale and is fitted with
[`stats::lm()`](https://rdrr.io/r/stats/lm.html), and the limits from
[`predict()`](https://rdrr.io/r/stats/predict.html) carry no uncertainty
in \\c\\. Fixing \\c\\ at the value estimated by a first fit, as in the
examples, treats an estimated \\c\\ as known: the curve is the same, but
the limits are narrower than they should be, because the uncertainty in
\\c\\ is left out.

## Variance

There is no `variance` argument. Equal scatter on the log scale already
means a standard deviation proportional to the flow, which is usually
why one would model the variance of a fit on the original scale, so it
is not implemented here. Taking logs does not always even out the
scatter completely, though, and a variance scheme on the log scale could
still be useful; check the residuals of the fit.

## Examples

``` r
fit <- rc_powerlaw_log(discharge, stage, data = thompson)
coef(fit)
#>         a         b         c 
#> 52.398628  1.880339 -1.304109 
predict(fit, new_stage = c(1, 3, 6), conflev = 0.95)
#> # A tibble: 3 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  252.   248.   255.
#> 2     3  815.   805.   825.
#> 3     6 2204.  2177.  2231.

# the offset known, say from a survey of the control
rc_powerlaw_log(discharge, stage, data = thompson, offset = -1.3)
#> Rating curve model.
#> - Method: rc_powerlaw_log 

# c estimated, then treated as known: the same curve, with narrower limits
# that leave out the uncertainty in c
fixed <- rc_powerlaw_log(
  discharge,
  stage,
  data = thompson,
  offset = fit$curve_parameters$c
)
predict(fixed, new_stage = c(1, 3, 6), conflev = 0.95)
#> # A tibble: 3 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  252.   248.   255.
#> 2     3  815.   808.   822.
#> 3     6 2204.  2177.  2230.
```
