# Contributing to CSHShydRometry

## How to contribute

Anyone interested is encouraged to contribute to the repository by
**forking and submitting a pull request**.

(If you are new to GitHub, you might start with a [basic
tutorial](https://docs.github.com/en/get-started/getting-started-with-git/set-up-git)
and check out a more detailed guide to [pull
requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).)

### Bug reports and feature requests

- Submit an issue on the [Issues
  page](https://github.com/CSHS-CWRA/CSHShydRometry/issues). For a bug,
  please include a small reproducible example and the output of
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html); or
- Submit a clean pull request with your code contribution for review.
  Pull requests will be reviewed by the maintainers and, if deemed
  beneficial, merged into `main`.

### Code contributions

- Fork this repo to your GitHub account.
- Clone your fork to your machine, e.g.
  `git clone https://github.com/[your username]/CSHShydRometry.git`
- Track progress upstream by adding this repository as a remote:
  `git remote add upstream https://github.com/CSHS-CWRA/CSHShydRometry.git`.
  Before making changes, pull in the latest from upstream, either with
  `git fetch upstream` and a later merge, or with
  `git pull upstream main`.
- Make your changes, preferably on a new branch.
- Before submitting, please:
  - add or update tests in `tests/testthat/` for any change in
    behaviour;
  - regenerate the documentation with `devtools::document()` if you
    changed any roxygen comments;
  - rebuild `README.md` with `devtools::build_readme()` if you changed
    `README.Rmd` (edit `README.Rmd`, not `README.md`);
  - add a line to `NEWS.md` describing any user-facing change;
  - check that `devtools::check()` passes.
- Push to your fork, and submit a pull request to `main` at
  `CSHS-CWRA/CSHShydRometry`.

### Roadmap

Directions the package could grow in. Each is open to contributors: if
you would like to take one on, open an issue first, so that the design
can be agreed before the work starts. Items marked *(to discuss)* are
undecided.
[DESIGN.md](https://cshs-cwra.github.io/CSHShydRometry/DESIGN.md) has
the reasoning behind several of them.

**Models**

- Curves with more than two segments, as `rc_pseg_*()`. The search over
  starting breakpoints would need to cover several breakpoints, kept in
  order.
- More two-segment fits: first, on the log-log scale
  (`rc_2seg_powerlaw_log()`). *(to discuss)* Fits whose segments are not
  power laws; hydraulic controls are usually modelled as power laws, so
  this may not be needed.
- A known offset for each segment in
  [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md),
  as
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md)
  and
  [`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md)
  have.
- More ways for segments to combine than `combine = "replace"` and
  `"add"`, which multi-segment curves will need. BaRatin describes them
  with a configuration matrix: a row for each range of stage, a column
  for each segment, and a 1 where that segment carries flow in that
  range, so “replace” is `rbind(c(1, 0), c(0, 1))` and “add” is
  `rbind(c(1, 0), c(1, 1))`. `combine` could accept such a matrix,
  keeping “replace” and “add” as shorthands, but users should never have
  to write one: common configurations need names.
  [`rc_2seg_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_2seg_powerlaw.md)
  stays as a convenience for two segments. *(to discuss)* How to name
  the configurations; see the design notes.
- A breakpoint-averaged interval method for two-segment curves, which
  removes the jump in the delta-method band at the breakpoint.
- *(to discuss)* Bayesian fitting.

**Variance**

- Estimating the power in
  [`var_power()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  within the package, for example by maximising the profile likelihood
  over the power, refitting at each value by reweighting in rounds. That
  would make it available in every fitting function, not only
  [`rc_powerlaw()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw.md),
  and remove the dependence on
  [`nlme::gnls()`](https://rdrr.io/pkg/nlme/man/gnls.html).
- Further estimated shapes for the scatter, such as an exponential
  (`var_exp()`).
- *(priority)* An error model that combines each gauging’s reported
  uncertainty with extra scatter estimated from the fit, as BaRatin does
  (its “remnant error”). This is the usual error model in the field, and
  would give prediction limits where
  [`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
  cannot. Variance schemes could be combined by adding them, such as
  `var_prop() + var_spec()`.
- Variance schemes on the log-log scale, in
  [`rc_powerlaw_log()`](https://cshs-cwra.github.io/CSHShydRometry/reference/rc_powerlaw_log.md).

**Predicting and uncertainty**

- A bootstrap for single-segment curves, and with it a `method` argument
  for their [`predict()`](https://rdrr.io/r/stats/predict.html) methods,
  as two-segment fits have. Until then they have no `method` argument;
  passing one is an error, as is any argument a method does not take.
- Bootstrap settings (`B`, `seed`, …) given as an object, like the
  variance schemes, rather than passed through `...`.
- Predicting a distribution: the predictive distribution of discharge at
  each stage, as a probaverse distribution, rather than only limits. The
  limits [`predict()`](https://rdrr.io/r/stats/predict.html) returns
  would be quantiles of it.
- *(to discuss)* The default interval method for multi-segment curves.
- A faster bootstrap, running resamples in parallel.

**Usability**

- A [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method,
  drawing the curve with stage on the vertical axis.
- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) methods.
- *(priority)* Ratings that shift over time, as the Thompson gaugings
  do: checking for drift, fitting by period, and curves that change with
  time. This would add an optional time argument (such as `date`) to the
  fitting functions and a matching `new_date` to
  [`predict()`](https://rdrr.io/r/stats/predict.html), with a `date`
  column in its output.

**Data**

- A dataset with a clear change of control, under a licence that allows
  it to ship with the package (the Sauze gaugings come from RBaM, which
  is GPL-3).

### Design decisions

Before changing how something works, read
[DESIGN.md](https://cshs-cwra.github.io/CSHShydRometry/DESIGN.md): it
records the decisions behind the package’s design, why they were made,
and what they leave open. If your change revisits one of them, update
its entry.

### Code style

- Match the style of the surrounding code.
- Name things in full rather than with terse abbreviations: `discharge`
  and `stage`, not `q` and `h`.
- In exported functions, put `...` straight after the required
  arguments, so that every optional argument has to be named in full.
- `discharge` and `stage` are looked up in `data`, so users can refer to
  its columns directly. Arguments of the `var_*()` functions are
  ordinary arguments; see
  [DESIGN.md](https://cshs-cwra.github.io/CSHShydRometry/DESIGN.md).

### Prefer to email?

Email the person listed as maintainer in the `DESCRIPTION` file of this
repo.

### Thanks for contributing!

Please note that this project is released with a [Contributor Code of
Conduct](https://cshs-cwra.github.io/CSHShydRometry/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.
