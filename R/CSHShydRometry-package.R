#' @keywords internal
"_PACKAGE"

#' @section Fitting a rating curve:
#' Every model is fitted by a `rc_*()` constructor and summarised by a
#' `predict()` method. The constructors take discharge and stage — as vectors,
#' or as column names with a `data` argument — and return an object carrying
#' `"rating_curve"` as its second class.
#'
#' Single segment:
#' \itemize{
#'   \item [rc_log_ols()], [rc_log_nls()] — power law fitted on the log-log
#'     scale, by ordinary or nonlinear least squares.
#'   \item [rc_nls()] — power law fitted on the natural scale.
#'   \item [rc_gnls()] — power law by generalised nonlinear least squares,
#'     with the error variance estimated as a power of the mean.
#'   \item [rc_poly()], [rc_loess()] — polynomial and loess alternatives.
#' }
#'
#' Two segments, joined at an estimated breakpoint:
#' \itemize{
#'   \item [rc_nls_2seg()] — `config = "piecewise"` forces the two power laws
#'     to meet at the breakpoint; `config = "compound"` adds the upper segment
#'     to the discharge carried at the breakpoint.
#' }
#'
#' @section Weighting:
#' The constructors that take `wts_code` offer three error models: `"none"`
#' (constant variance), `"prop"` (constant coefficient of variation, fitted by
#' iteratively reweighted least squares), and `"spec"` (variances supplied by
#' the user, typically from reported gauging uncertainties). Under `"spec"` a
#' new observation's scatter is not identified by the fit, so prediction limits
#' are returned as `NA`.
#'
#' @section Confidence and prediction limits:
#' `predict()` returns the same columns whatever the model, the method or the
#' weighting: `h`, `fit`, and — when `conflev` or `predlev` is given —
#' `ci_lwr`/`ci_upr` and `pi_lwr`/`pi_upr`. Quantities that cannot be computed
#' come back as `NA` rather than as missing columns, so results from different
#' approaches stack directly with `rbind()`.
#'
#' For two-segment curves, [predict.rc_nls_2seg()] takes a `method`. The
#' default, `"delta"`, is fast but unreliable near the breakpoint, where the
#' mean function is not differentiable; `"boot"` costs a refit per resample and
#' behaves much better there. See [predict.rc_nls_2seg()] for the detail.
NULL
