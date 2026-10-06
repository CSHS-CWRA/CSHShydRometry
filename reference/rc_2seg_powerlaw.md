# Fit two-segment power-law rating curve using nls on untransformed data

Fit two-segment power-law rating curve using nls on untransformed data

## Usage

``` r
rc_2seg_powerlaw(
  discharge,
  stage,
  ...,
  data = NULL,
  combine = c("replace", "add"),
  kstart = NULL,
  kfixed = FALSE,
  kbounds = NULL,
  variance = var_none(),
  control = stats::nls.control(maxiter = 1000)
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

  Must be empty. Present so that every argument after it has to be
  named, which keeps calls readable and guards against positional
  mistakes.

- data:

  Optional data frame in which `discharge` and `stage` are evaluated.

- combine:

  How the two segments combine above the breakpoint. `"replace"` (the
  default): the upper power law takes over from the lower one, with its
  coefficient set so that the two meet at the breakpoint. `"add"`: the
  upper power law adds to the discharge the lower one carries at the
  breakpoint, as when flow spills onto a floodplain.

- kstart:

  Starting value(s) for the breakpoint `k`. `NULL`, the default, tries
  10 values spread evenly across the search range; a numeric vector
  tries each of its values. See Details.

- kfixed:

  If `TRUE`, hold the breakpoint `k` fixed at `kstart`.

- kbounds:

  Lower and upper bounds for `k`, or `NULL` to keep at least three
  gaugings in each segment.

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

- control:

  Settings for [`stats::nls()`](https://rdrr.io/r/stats/nls.html), as
  from
  [`stats::nls.control()`](https://rdrr.io/r/stats/nls.control.html).
  The fit uses the `"port"` algorithm, which ignores `tol`; its own
  settings, such as `rel.tol`, can be added to the list (see
  [`stats::nls()`](https://rdrr.io/r/stats/nls.html)).

## Value

An object of class `c("rc_2seg_powerlaw", "rating_curve")`; see
[rating_curve](https://cshs-cwra.github.io/CSHShydRometry/reference/rating_curve.md)
for its contents. `curve_parameters` holds the estimated parameters by
type, one value per segment: under `"replace"`, `a` has a single value
because the upper segment's coefficient is fixed by continuity, and
under `"add"`, `c` has a single value because the upper segment is
measured from `k`. [`coef()`](https://rdrr.io/r/stats/coef.html) gives
the same estimates by their model names (`a1`, `b1`, `c1`, ...).
`kstart_search` is a tibble with a row per starting breakpoint tried:
the estimated breakpoint `k` it led to and the `loss` of that fit (see
Details), both `NA` where the fit failed.

## Details

Fitting proceeds in three steps, the last two repeated for each starting
breakpoint:

1.  Choose the breakpoint search range (from `kbounds`, or a default
    that leaves at least three gaugings in each segment) and the
    starting breakpoints to try.

2.  Derive starting values for the segment parameters by splitting the
    data at the starting breakpoint and fitting a linear model to each
    segment on the log-log scale (`log(discharge) ~ log(stage - c)`);
    `a = exp(intercept)`, `b = slope`.

3.  Fit all parameters jointly with
    [`stats::nls()`](https://rdrr.io/r/stats/nls.html) using the "port"
    algorithm (which supports the box constraints in `lower`/`upper`).
    Under
    [`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
    this fit is repeated, updating the weights from the current fitted
    values and starting from the previous estimates, until the fitted
    discharges stabilise.

The two-segment fit is sensitive to where it starts: from some starting
breakpoints [`nls()`](https://rdrr.io/r/stats/nls.html) fails outright,
and from others it settles on a local optimum. Trying several starts
guards against both. Of the fits that succeed, the one kept has the
smallest loss, the quantity the fit minimises. Under
[`var_none()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
and
[`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
that is the weighted residual sum of squares. Under
[`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
it is not what you might expect: the reweighting does not minimise the
sum of squared relative residuals, but settles where the Gamma
quasi-likelihood is maximised, so the loss is
`sum(discharge / fitted + log(fitted))`. See Wedderburn (1974) for
quasi-likelihood. Under
[`var_prop()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md),
fits whose reweighting converged are preferred. What each start led to
is recorded in `kstart_search`.

Each start costs a full fit, and the bootstrap
(`predict(method = "boot")`) repeats the same search for every resample.
That is deliberate: a resample's best fit is often found from a
different start than the original's, so starting the resamples only from
the original estimate, or only from the optima the original search
found, would understate the uncertainty in the breakpoint. A single
`kstart` makes both the fit and the bootstrap faster, but both then rest
on that one start.

## References

Wedderburn, R. W. M. (1974). Quasi-likelihood functions, generalized
linear models, and the Gauss-Newton method. *Biometrika*, 61(3),
439-447.
[doi:10.1093/biomet/61.3.439](https://doi.org/10.1093/biomet/61.3.439)

## Examples

``` r
# The Thompson is close to a single control, so its two-segment fit needs
# proportional weights to converge.
fit <- rc_2seg_powerlaw(discharge, stage, data = thompson, variance = "prop")
fit
#> Rating curve model.
#> - Method: rc_2seg_powerlaw 
coef(fit)
#>          a1          b1          c1          b2          c2           k 
#> 192.5532263   0.9018095  -0.3288885   1.7185012  -0.8851584   1.1800577 

# what each starting breakpoint led to
fit$kstart_search
#> # A tibble: 10 × 3
#>    kstart     k  loss
#>     <dbl> <dbl> <dbl>
#>  1   1.19 NA      NA 
#>  2   1.85  1.18  716.
#>  3   2.50 NA      NA 
#>  4   3.16  1.18  716.
#>  5   3.82 NA      NA 
#>  6   4.48 NA      NA 
#>  7   5.14 NA      NA 
#>  8   5.80 NA      NA 
#>  9   6.46 NA      NA 
#> 10   7.12 NA      NA 

predict(fit, new_stage = c(1, 3, 6), conflev = 0.95)
#> # A tibble: 3 × 4
#>   stage   fit ci_lwr ci_upr
#>   <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  249.   243.   255.
#> 2     3  827.   814.   839.
#> 3     6 2210.  2184.  2236.

# A river with a clearer change of control is far less fussy. The Ardeche
# at Sauze, in the RBaM package, fits however the segments combine, and
# carries a reported uncertainty for every gauging.
if (requireNamespace("RBaM", quietly = TRUE)) {
  sauze <- RBaM::SauzeGaugings
  repl <- rc_2seg_powerlaw(Q, H, data = sauze, kstart = 1)
  add <- rc_2seg_powerlaw(Q, H, data = sauze, combine = "add", kstart = 1)
  c(replace = repl$curve_parameters[["k"]], add = add$curve_parameters[["k"]])

  # variances from the reported gauging uncertainties
  rc_2seg_powerlaw(
    Q,
    H,
    data = sauze,
    variance = var_spec(sauze$uQ^2),
    kstart = 1
  )
}
#> Rating curve model.
#> - Method: rc_2seg_powerlaw 
```
