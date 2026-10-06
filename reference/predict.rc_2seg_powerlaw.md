# Predict discharge with confidence / prediction limits from a two-segment fit

Evaluates the fitted rating curve at the stages given and, optionally,
attaches confidence and/or prediction limits, by the delta method or a
bootstrap.

## Usage

``` r
# S3 method for class 'rc_2seg_powerlaw'
predict(
  object,
  ...,
  new_stage = NULL,
  conflev = NULL,
  predlev = NULL,
  method = c("delta", "boot")
)
```

## Arguments

- object:

  An `rc_2seg_powerlaw` fit (from
  [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md)).

- ...:

  For `method = "boot"`, its settings:

  `B`

  :   The number of resamples that must fit successfully (default 1000).

  `seed`

  :   An optional random seed, for reproducible limits.

  `max_tries_factor`

  :   A cap on the attempts, `B * max_tries_factor` (default 3), so that
      the bootstrap ends even if many resamples fail.

  For `method = "delta"`, must be empty. Any other argument is an error.

- new_stage:

  Stages at which to return limits. Defaults to 1000 points spanning the
  observed stage range.

- conflev:

  Confidence level for the mean-curve (confidence) interval, or `NULL`
  to omit it.

- predlev:

  Confidence level for the prediction interval, or `NULL` to omit it.

- method:

  Interval method:

  - `"delta"` (the default): the linearised delta method. Fast, but
    fixes the breakpoint at its estimate and so jumps there; see
    "Choosing a method".

  - `"boot"`: a bootstrap, which resamples the gaugings and refits; see
    "The bootstrap". Much the slowest, since it refits the model `B`
    times, and the most trustworthy at the breakpoint.

## Value

A tibble with column `stage`, the fitted discharge `fit`, and, when
requested, `ci_lwr`/`ci_upr` and `pi_lwr`/`pi_upr`. Which columns are
present depends only on which of `conflev` and `predlev` were given –
never on the method or the weighting. Where a quantity cannot be
computed (prediction limits under `"spec"` weights) the column is
returned as `NA`.

Named for the object's class, `rc_2seg_powerlaw`, so `predict(object)`
dispatches here. The one-segment models define their own
`predict.rc_powerlaw` for the `rc_powerlaw` class; the two do not
collide.

## Choosing a method

The default, `"delta"`, is fast and is the classical choice, but **it is
not reliable near the breakpoint**. The mean function is not
differentiable at `stage = k`, so the linearisation switches form there
and the interval jumps: on one fitted curve the band widened from 29 to
242 m^3 s^-1 across the breakpoint. In a simulation study a nominal 95%
delta interval covered the true curve only about two-thirds of the time
just above the breakpoint, against roughly 97% for `"boot"`. Away from
the breakpoint it behaves normally.

So: `"delta"` for a quick look or where the breakpoint is not of
interest, and `"boot"` where the interval matters.

A third approach – drawing parameter vectors from their asymptotic
normal distribution and pushing each draw through the model – was tried
and dropped. When a segment is poorly identified, as the lower segment
at Sauze is, the draws too often land on impossible curves: negative or
astronomically large discharges, giving limits that are meaningless.

## The bootstrap

`method = "boot"` resamples the gaugings with replacement, refits the
two-segment curve to each resample with the arguments the fit was made
with (repeating the search over starting breakpoints), and summarises
the spread of the refitted curves. Under
[`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
the variances are resampled along with the gaugings. A resample fails if
its fit errors or, under
[`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
if its reweighting does not converge; failed resamples are redrawn, up
to `B * max_tries_factor` attempts in all.

Confidence limits are the pointwise percentiles of the refitted curves.
Prediction limits add the scatter of a new gauging to the spread of the
refitted curves, using the t distribution on the fit's residual degrees
of freedom. Under
[`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
they are `NA`, as for the delta method. The number of resamples that
fitted is recorded in `attr(, "B_success")`.

## What the limits assume

Prediction limits, from either method, assume the scatter of the
gaugings about the curve is normal, with the spread the weighting scheme
describes: if it is skewed or has heavy tails, they can miss, especially
at high levels such as 0.99. Delta-method confidence limits assume the
estimates are close to normal, which is approximate, and poor near the
breakpoint; bootstrap confidence limits do not.

## Examples

``` r
if (requireNamespace("RBaM", quietly = TRUE)) {
  sauze <- RBaM::SauzeGaugings
  fit <- rc_2seg_powerlaw(Q, H, data = sauze, kstart = 1)
  hp <- c(1, 1.5, 2, 4)

  # the default, and fast
  predict(fit, new_stage = hp, conflev = 0.95)

  # slower, but does not assume the breakpoint is known. Compare the two
  # either side of the breakpoint, at about 1.85 m: the delta band jumps
  # there, the bootstrap band does not.
  predict(fit, new_stage = hp, conflev = 0.95, method = "boot", B = 50)
}
#> # A tibble: 4 × 4
#>   stage    fit ci_lwr ci_upr
#>   <dbl>  <dbl>  <dbl>  <dbl>
#> 1   1     84.0   78.5   95.3
#> 2   1.5  141.   120.   184. 
#> 3   2    259.   201.   308. 
#> 4   4   1107.  1004.  1290. 

# the columns returned never depend on the model or the method, so results
# from different approaches stack directly
one <- rc_powerlaw(discharge, stage, data = thompson)
two <- rc_2seg_powerlaw(discharge, stage, data = thompson, variance = "prop", kstart = 2)
rbind(
  predict(one, new_stage = 3, conflev = 0.95),
  predict(two, new_stage = 3, conflev = 0.95)
)
#> # A tibble: 2 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     3  832.   814.   850.
#> 2     3  827.   814.   839.
```
