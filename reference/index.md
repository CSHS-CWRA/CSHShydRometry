# Package index

## Fitting a single-segment curve

Each function fits one model, and returns a fit that
[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`coef()`](https://rdrr.io/r/stats/coef.html) work with.

- [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md)
  : Fit a power-law rating curve on the original scale
- [`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md)
  : Fit a power-law rating curve on the log-log scale
- [`rc_poly()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_poly.md)
  : Fit polynomial rating curve
- [`rc_loess()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_loess.md)
  : Fit rating curve using loess smoother

## Fitting a two-segment curve

- [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md)
  : Fit two-segment power-law rating curve using nls on untransformed
  data

## Variance

How the scatter of the gaugings about the curve is modelled.

- [`var_none()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  [`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  [`print(`*`<rc_variance>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  : Variance schemes for rating-curve fits

## Predicting, with confidence and prediction limits

- [`predict(`*`<rc_2seg_powerlaw>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/predict.rc_2seg_powerlaw.md)
  : Predict discharge with confidence / prediction limits from a
  two-segment fit
- [`predict(`*`<rc_loess>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/predict.rc_loess.md)
  : Predict method for rc_loess objects
- [`predict(`*`<rc_poly>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/predict.rc_poly.md)
  : Predict method for rc_poly objects
- [`predict(`*`<rc_powerlaw>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/predict.rc_powerlaw.md)
  : Predict method for rc_powerlaw objects
- [`predict(`*`<rc_powerlaw_log>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/predict.rc_powerlaw_log.md)
  : Predict method for rc_powerlaw_log objects

## Fitted curves

What a fit contains, and the methods every fit shares.

- [`rating_curve`](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md)
  : Rating curve fits
- [`coef(`*`<rating_curve>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/coef.rating_curve.md)
  [`coef(`*`<rc_2seg_powerlaw>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/coef.rating_curve.md)
  : Estimated parameters of a rating curve
- [`fitted(`*`<rating_curve>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/fitted.rating_curve.md)
  [`residuals(`*`<rating_curve>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/fitted.rating_curve.md)
  : Fitted values and residuals of a rating curve
- [`print(`*`<rating_curve>`*`)`](https://cshs-cwra.github.io/CSHShydRometry/reference/print.rating_curve.md)
  : Print a rating curve fit

## Data

- [`thompson`](https://cshs-cwra.github.io/CSHShydRometry/reference/thompson.md)
  : Thompson River gaugings
