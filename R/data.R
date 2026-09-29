# Documentation for the datasets shipped with the package.

#' Thompson River gaugings
#'
#' Stage and discharge measurements from Water Survey of Canada station
#' 08LF051, Thompson River, used in the examples. The raw export carries one
#' row per field activity, including stage-only visits; this is the subset with
#' a discharge measurement, reduced to the columns a rating curve needs and
#' sorted by stage.
#'
#' Reported uncertainties are available for only a minority of the gaugings, so
#' `uq` is mostly `NA`. Fits using `wts_code = "spec"` need a complete weight
#' vector, and so need either a subset with `uq` present or weights from
#' elsewhere.
#'
#' @format A data frame with 93 rows and 4 columns:
#' \describe{
#'   \item{date}{Date of the measurement.}
#'   \item{h}{Stage, in metres (the mean gauge height for the activity).}
#'   \item{q}{Discharge, in cubic metres per second.}
#'   \item{uq}{Reported discharge uncertainty, or `NA` where none was given.}
#' }
#' @source Water Survey of Canada, station 08LF051. Built by
#'   `data-raw/thompson.R`.
"thompson"
