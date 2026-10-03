# The package-level help page, ?CSHShydRometry. Its title, description and
# links come from DESCRIPTION; the sections below are added to them. The
# authors are written out because roxygen's default lists the maintainer
# separately from, and not among, the authors: keep them in step with
# Authors@R in DESCRIPTION.

#' @section Fitting a rating curve:
#' Every model is fitted by a `rc_*()` constructor and evaluated by a
#' `predict()` method. The constructors take discharge and stage — as vectors,
#' or as column names with a `data` argument — and return an object carrying
#' `"rating_curve"` as its top-level class; see [rating_curve].
#'
#' Single segment:
#' \itemize{
#'   \item [rc_powerlaw()] — power law fitted by least squares on the original
#'     scale.
#'   \item [rc_powerlaw_log()] — power law fitted by least squares on the
#'     log-log scale.
#' }
#' In both, the offset \eqn{c} is estimated, or held at a value given as
#' `offset`.
#' \itemize{
#'   \item [rc_poly()], [rc_loess()] — polynomial and loess alternatives.
#' }
#'
#' Two segments, joined at an estimated breakpoint:
#' \itemize{
#'   \item [rc_2seg_powerlaw()] — with `combine = "replace"`: the upper power
#'     law takes over from the lower at the breakpoint; with 
#'     `combine = "add"`: it adds to the discharge carried at the breakpoint.
#' }
#'
#' @section Weighting:
#' The constructors that take `wts` offer four models of how the gaugings
#' scatter about the curve:
#' \itemize{
#'   \item [wts_none()] — the same scatter at every flow (constant variance).
#'   \item [wts_prop()] — scatter proportional to the flow.
#'   \item [wts_power()] — scatter proportional to a power of the flow, with
#'     the power estimated; [rc_powerlaw()] only at this time.
#'   \item [wts_spec()] — variances supplied by the user, typically from
#'     reported gauging uncertainties. A new observation's scatter is then
#'     not identified by the fit, so prediction limits are returned as `NA`.
#' }
#' See [wts].
#'
#' @section Confidence and prediction limits:
#' `predict()` returns the same columns whatever the model, the method or the
#' weighting: `stage`, `fit`, and — when `conflev` or `predlev` is given —
#' `ci_lwr`/`ci_upr` and `pi_lwr`/`pi_upr`. Quantities that cannot be computed
#' come back as `NA` rather than as missing columns, so results from different
#' approaches stack directly.
#'
#' For two-segment curves, [predict.rc_2seg_powerlaw()] takes a `method`. The
#' default, `"delta"`, is fast but unreliable near the breakpoint, where the
#' mean function is not differentiable; `"boot"` costs a refit per resample
#' and behaves much better there. See [predict.rc_2seg_powerlaw()] for the detail.
#'
#' @section Data:
#' [thompson] ships with the package and suits the single-segment models. The
#' two-segment examples prefer the Ardeche at Sauze, `RBaM::SauzeGaugings`,
#' which has a clearer change of control; RBaM is a suggested dependency used
#' only as a source of that data.
#'
#' @author
#' Vincenzo Coia (maintainer,
#' \email{vincenzo.coia@@gmail.com}), R. Dan Moore, and Paul Whitfield.
#' @keywords internal
"_PACKAGE"

# Import from tibble so that loading this package loads it. Every table the
# package returns or ships is a tibble, and a tibble only behaves like one
# (printing, `[`, no partial matching with `$`) once tibble is loaded;
# calling it with `tibble::` alone would load it only on first use.
#' @importFrom tibble tibble
NULL
