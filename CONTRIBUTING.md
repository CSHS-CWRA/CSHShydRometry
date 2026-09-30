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
