# Rating curve fits

Every `rc_*()` constructor returns a list with class
`c("rc_<method>", "rating_curve")`. These elements are part of the
interface:

## Details

- `gaugings`:

  A tibble of the gaugings used, with columns `discharge` and `stage`,
  after dropping any with a missing value.

- `curve_parameters`:

  The estimated parameters of the curve, as a named list with one
  element per parameter type. For a multi-segment curve an element holds
  one value per segment (for `k`, per breakpoint); elements may differ
  in length where a parameter is fixed by how the segments join. Empty
  for
  [`rc_loess()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_loess.md),
  which has no parameters. [`coef()`](https://rdrr.io/r/stats/coef.html)
  returns the same estimates as a flat named vector.

- `settings`:

  The arguments the fit was made with, including any user-supplied
  variances, aligned with `gaugings`. Enough to refit.

- `variance`:

  For constructors with a `variance` argument: the variance scheme, with
  any estimated parameter filled in, such as the power estimated under
  [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  (`fit$variance$exponent`).

Use [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) rather than reaching
into the fit for those quantities.

Fits also carry elements that depend on how they are computed, and that
may change as the methods for modelling the scatter develop: `model`,
the underlying model object, such as from
[`stats::nls()`](https://rdrr.io/r/stats/nls.html); `rse`, the residual
standard error (on the log scale for the log-scale fits); and, for
constructors with a `variance` argument, `weights_used`, the weights the
final model was fitted with, and `irls` (for iteratively reweighted
least squares): `NULL`, or under
[`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
which refits in rounds, a list giving the number of rounds
(`iterations`) and whether they `converged`. Some constructors carry
further elements, documented on their own help pages.
