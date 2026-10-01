# Build the `thompson` dataset shipped with the package.
#
# Source: Water Survey of Canada gauging record for station 08LF051,
# Thompson River near Spences Bridge, as received from the Water Survey
# (field-visit records are not published on its Water Level and Flow
# website). The raw export carries one
# row per field activity, including stage-only visits and a large number of
# administrative columns; only the discharge measurements are kept, and only
# the columns a rating curve needs.
#
#   Rscript data-raw/thompson.R

raw_csv <- file.path("data-raw", "08LF051_rating_data.csv")
if (!file.exists(raw_csv)) {
  stop("missing ", raw_csv)
}

raw <- utils::read.csv(raw_csv, check.names = TRUE, stringsAsFactors = FALSE)

keep <- !is.na(raw$Discharge) &
  !is.na(raw$Mean.Gauge.Height) &
  raw$Discharge > 0

# The reported uncertainty is a percentage of the discharge, at two standard
# deviations: the activity remarks describe it as the "IVE method, 2-sigma
# value". It is kept as reported.
thompson <- data.frame(
  date = as.Date(substr(raw$Date..UTC.[keep], 1, 10)),
  stage = as.numeric(raw$Mean.Gauge.Height[keep]),
  discharge = as.numeric(raw$Discharge[keep]),
  uncertainty_pct = as.numeric(raw$Uncertainty[keep])
)
thompson <- thompson[order(thompson$stage), ]
thompson <- tibble::as_tibble(thompson)

stopifnot(
  nrow(thompson) > 50,
  all(thompson$discharge > 0),
  !any(is.na(thompson$stage)),
  !any(is.na(thompson$discharge))
)

save(thompson, file = file.path("data", "thompson.rda"), compress = "bzip2")
cat("wrote data/thompson.rda:", nrow(thompson), "gaugings\n")
