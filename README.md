# CSHShydRometry

Frequentist methods for fitting stage–discharge rating curves, with confidence
and prediction limits.

Derived from work by Dan Moore on rating-curve methods, restructured into a
model–predict form: every model is fitted by a `rc_*()` constructor and
summarised by a `predict()` method.

## Installation

``` r
# install.packages("pak")
pak::pak("CSHShydRometry")
```

## Fitting a curve

``` r
library(CSHShydRometry)

fit <- rc_nls(q, h, data = thompson)
fit
#> Rating curve model.
#> - Method: rc_nls

predict(fit, hpred = c(1, 3, 6), conflev = 0.95)
```

Single-segment models: `rc_log_ols()` and `rc_log_nls()` fit a power law on the
log–log scale, `rc_nls()` on the natural scale, `rc_gnls()` estimates the error
variance as a power of the mean, and `rc_poly()` and `rc_loess()` offer
non-power-law alternatives.

Two-segment curves are fitted by `rc_nls_2seg()`, joined either continuously
(`config = "piecewise"`) or additively (`config = "compound"`) at an estimated
breakpoint.

## Consistent output

`predict()` returns the same columns whatever the model, method or weighting:
`h`, `fit`, and — when asked for — `ci_lwr`/`ci_upr` and `pi_lwr`/`pi_upr`.
Quantities that cannot be computed come back as `NA` rather than as missing
columns, so a batch of approaches stacks directly:

``` r
rbind(
  predict(rc_nls(q, h, data = thompson),  hpred = 3, conflev = 0.95),
  predict(rc_poly(q, h, data = thompson), hpred = 3, conflev = 0.95),
  predict(rc_loess(q, h, data = thompson), hpred = 3, conflev = 0.95)
)
```

## Weighting

Constructors taking `wts_code` offer three error models:

| `wts_code` | assumption |
|---|---|
| `"none"` | constant variance |
| `"prop"` | constant coefficient of variation, fitted by IRLS |
| `"spec"` | variances supplied by the user, e.g. from reported gauging uncertainties |

Under `"spec"` a new observation's scatter is not identified by the fit, so
prediction limits are returned as `NA`.

## Intervals for two-segment curves

`predict.rc_nls_2seg()` takes a `method`:

| method | notes |
|---|---|
| `"delta"` (default) | linearised, fast, **unreliable near the breakpoint** |
| `"boot"` | resamples the gaugings and refits; slow, but trustworthy at the breakpoint |
| `"sim"` | draws parameters from their asymptotic normal distribution; sensitive to the parameter scale |

The default warrants a word. A two-segment mean function is not differentiable
at the breakpoint, so the delta method's linearisation switches form there and
the interval jumps. On one fitted curve the band widened from 29 to
242 m³ s⁻¹ across the breakpoint, and in a simulation study a nominal 95%
delta interval covered the true curve only about two-thirds of the time just
above it, against roughly 97% for the bootstrap. Away from the breakpoint the
delta method behaves normally.

So `"delta"` for a quick look, and `"boot"` where the interval matters.

## Scope

This package covers the frequentist methods only. Bayesian rating-curve
estimation, and a breakpoint-averaged interval method that removes the jump
described above, are deliberately left out for now.
