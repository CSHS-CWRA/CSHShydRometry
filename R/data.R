# Documentation for the datasets shipped with the package.

#' Thompson River gaugings
#'
#' Stage and discharge measurements from Water Survey of Canada station
#' 08LF051, Thompson River, used in the examples: every discharge measurement
#' made at the station from 1994 to 2024, sorted by stage.
#'
#' Reported uncertainties are available for only 19 of the gaugings, so
#' `uncertainty_pct` is mostly `NA`. It is a percentage of the discharge at
#' two standard deviations, so the standard uncertainty of a discharge, in
#' cubic meters per second, is `uncertainty_pct / 100 * discharge / 2`.
#' Fits using `wts_code = "spec"` need a complete weight vector, and so need
#' either the subset with an uncertainty or weights from elsewhere.
#'
#' One uncertainty looks wrong: the gauging of 2018-11-01 has 0.0226, where
#' the others lie between about 2.5 and 11. It was probably entered as a
#' fraction (2.26%) rather than a percentage. The value is like this in the
#' file downloaded from the Water Survey, so the slip happened at the source,
#' not in this package. It is kept as reported.
#'
#' @section Getting the data yourself:
#' The measurements are published by the Water Survey of Canada on its Water
#' Level and Flow website, <https://wateroffice.ec.gc.ca>: search for station
#' 08LF051. The download there has more rows (visits that recorded only the
#' stage) and more columns than `thompson`.
#'
#' @format A data frame with 93 rows and 4 columns:
#' \describe{
#'   \item{date}{Date of the measurement.}
#'   \item{stage}{Stage, in meters (the mean gauge height for the visit).}
#'   \item{discharge}{Discharge, in cubic meters per second.}
#'   \item{uncertainty_pct}{Reported uncertainty of the discharge, as a
#'     percentage of it at two standard deviations (the "IVE method, 2-sigma
#'     value"), or `NA` where none was given.}
#' }
#' @source Water Survey of Canada, station 08LF051,
#'   <https://wateroffice.ec.gc.ca>.
"thompson"
