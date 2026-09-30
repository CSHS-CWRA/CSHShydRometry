# Documentation for the datasets shipped with the package.

#' Thompson River gaugings
#'
#' Stage and discharge measurements from Water Survey of Canada station
#' 08LF051, Thompson River, used in the examples. The raw export carries one
#' row per field activity, including stage-only visits; this is the subset with
#' a discharge measurement, reduced to the columns a rating curve needs and
#' sorted by stage.
#'
#' Reported uncertainties are available for only 19 of the gaugings, so
#' `uncertainty_pct` is mostly `NA`. It is a percentage of the discharge at
#' two standard deviations, so the standard uncertainty of a discharge, in
#' cubic metres per second, is `uncertainty_pct / 100 * discharge / 2`.
#' Fits using `wts_code = "spec"` need a complete weight vector, and so need
#' either the subset with an uncertainty or weights from elsewhere.
#'
#' One value, for the gauging of 2018-11-01, is 0.0226 where the others lie
#' between about 2.5 and 11; it may have been recorded as a fraction rather
#' than a percentage. It is in the export as downloaded, on the one record
#' whose only remark is "Sources:Import", which suggests it arrived from an
#' earlier system. It is kept as reported.
#'
#' @format A data frame with 93 rows and 4 columns:
#' \describe{
#'   \item{date}{Date of the measurement.}
#'   \item{stage}{Stage, in metres (the mean gauge height for the activity).}
#'   \item{discharge}{Discharge, in cubic metres per second.}
#'   \item{uncertainty_pct}{Reported uncertainty of the discharge, as a
#'     percentage of it at two standard deviations (the "IVE method, 2-sigma
#'     value"), or `NA` where none was given.}
#' }
#' @source Water Survey of Canada, station 08LF051. Built by
#'   `data-raw/thompson.R`.
"thompson"
