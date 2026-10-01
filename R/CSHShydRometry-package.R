#' @keywords internal
"_PACKAGE"

# Import from tibble so that loading this package loads it. Every table the
# package returns or ships is a tibble, and a tibble only behaves like one
# (printing, `[`, no partial matching with `$`) once tibble is loaded;
# calling it with `tibble::` alone would load it only on first use.
#' @importFrom tibble tibble
NULL

#' @section Fitting a rating curve:
#' Every model is fitted by a `rc_*()` constructor and evaluated by a
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
#'   \item [rc_2seg_nls()] — `controls = "successive"`: the upper power law
#'     takes over from the lower at the breakpoint; `controls = "additive"`: it adds
#'     to the discharge carried at the breakpoint.
#' }
#'
#' @section Weighting:
#' The constructors that take `wts` offer three error models: [wts_none()]
#' (constant variance), [wts_prop()] (constant coefficient of variation, fitted by
#' iteratively reweighted least squares), and [wts_spec()] (variances supplied
#' by the user, typically from reported gauging uncertainties). Under
#' [wts_spec()] a new observation's scatter is not identified by the fit, so
#' prediction limits are returned as `NA`. See [wts].
#'
#' @section Confidence and prediction limits:
#' `predict()` returns the same columns whatever the model, the method or the
#' weighting: `stage`, `fit`, and — when `conflev` or `predlev` is given —
#' `ci_lwr`/`ci_upr` and `pi_lwr`/`pi_upr`. Quantities that cannot be computed
#' come back as `NA` rather than as missing columns, so results from different
#' approaches stack directly with `rbind()`.
#'
#' @section Data:
#' [thompson] ships with the package and suits the single-segment models. The
#' two-segment examples prefer the Ardeche at Sauze, `RBaM::SauzeGaugings`,
#' which has a clearer change of control; RBaM is a suggested dependency used
#' only as a source of that data.
#'
#' For two-segment curves, [predict.rc_2seg_nls()] takes a `method`. The
#' default, `"delta"`, is fast but unreliable near the breakpoint, where the
#' mean function is not differentiable; `"boot"` costs a refit per resample and
#' behaves much better there. See [predict.rc_2seg_nls()] for the detail.
NULL
