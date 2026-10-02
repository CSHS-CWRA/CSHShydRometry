## How to contribute

Anyone interested is encouraged to contribute to the repository by **forking
and submitting a pull request**.

(If you are new to GitHub, you might start with a
[basic tutorial](https://docs.github.com/en/get-started/getting-started-with-git/set-up-git)
and check out a more detailed guide to
[pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests).)

### Bug reports and feature requests

* Submit an issue on the
  [Issues page](https://github.com/CSHS-CWRA/CSHShydRometry/issues). For a
  bug, please include a small reproducible example and the output of
  `sessionInfo()`; or
* Submit a clean pull request with your code contribution for review. Pull
  requests will be reviewed by the maintainers and, if deemed beneficial,
  merged into `main`.

### Code contributions

* Fork this repo to your GitHub account.
* Clone your fork to your machine, e.g.
  `git clone https://github.com/[your username]/CSHShydRometry.git`
* Track progress upstream by adding this repository as a remote:
  `git remote add upstream https://github.com/CSHS-CWRA/CSHShydRometry.git`.
  Before making changes, pull in the latest from upstream, either with
  `git fetch upstream` and a later merge, or with `git pull upstream main`.
* Make your changes, preferably on a new branch.
* Before submitting, please:
  * add or update tests in `tests/testthat/` for any change in behaviour;
  * regenerate the documentation with `devtools::document()` if you changed
    any roxygen comments;
  * rebuild `README.md` with `devtools::build_readme()` if you changed
    `README.Rmd` (edit `README.Rmd`, not `README.md`);
  * add a line to `NEWS.md` describing any user-facing change;
  * check that `devtools::check()` passes.
* Push to your fork, and submit a pull request to `main` at
  `CSHS-CWRA/CSHShydRometry`.

### Roadmap

Directions the package could grow in. Each is open to contributors: if you
would like to take one on, open an issue first, so that the design can be
agreed before the work starts. Items marked *(to discuss)* are undecided.
[DESIGN.md](DESIGN.md) has the reasoning behind several of them.

**Models**

* Curves with more than two segments, as `rc_pseg_*()`. The search over
  starting breakpoints would need to cover several breakpoints, kept in
  order.
* More two-segment fits: first, on the log-log scale
  (`rc_2seg_power_log()`). *(to discuss)* Fits whose segments are not power
  laws; hydraulic controls are usually modelled as power laws, so this may
  not be needed.
* A known stage of zero flow in `rc_power()` and `rc_2seg_power()`, as
  `rc_power_log()` has.
* A breakpoint-averaged interval method for two-segment curves, which
  removes the jump in the delta-method band at the breakpoint.
* *(to discuss)* Bayesian fitting.

**Weighting**

* Estimating the power in `wts_power()` within the package, for example by
  maximising the profile likelihood over the power, refitting at each value
  by reweighting in rounds. That would make it available in every fitting
  function, not only `rc_power()`, and remove the dependence on
  `nlme::gnls()`.
* Further estimated shapes, such as an exponential (`wts_exp()`), and a
  public `wts_nlme()` accepting any of nlme's variance functions.
* Known uncertainty plus extra scatter that grows with the flow, as one
  scheme (`wts_comb()`).
* Weighting on the log-log scale, in `rc_power_log()`.

**Predicting and uncertainty**

* `method` for every `predict()` method, not only two-segment fits, with a
  bootstrap for single-segment curves too.
* Bootstrap settings (`B`, `seed`, ...) given as an object, like the
  weighting schemes, rather than passed through `...`.
* *(to discuss)* The default interval method for multi-segment curves.
* A faster bootstrap, running resamples in parallel.

**Usability**

* A `plot()` method, drawing the curve with stage on the vertical axis.
* `fitted()` and `residuals()` methods.
* Ratings that shift over time: checking for drift, and fitting by period.

**Data**

* A dataset with a clear change of control, under a licence that allows it
  to ship with the package (the Sauze gaugings come from RBaM, which is
  GPL-3).

### Design decisions

Before changing how something works, read [DESIGN.md](DESIGN.md): it records
the decisions behind the package's design, why they were made, and what they
leave open. If your change revisits one of them, update its entry.

### Code style

* Match the style of the surrounding code.
* Name things in full rather than with terse abbreviations: `discharge` and
  `stage`, not `q` and `h`.
* In exported functions, put `...` straight after the required arguments, so
  that every optional argument has to be named in full.
* Per-gauging arguments (discharge, stage, weights) are looked up in `data`,
  so users can refer to its columns directly.

### Prefer to email?

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

### Thanks for contributing!

Please note that this project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this
project you agree to abide by its terms.
