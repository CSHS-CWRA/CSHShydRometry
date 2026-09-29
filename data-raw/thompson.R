# Build the `thompson` dataset shipped with the package.
#
# Source: Water Survey of Canada gauging record for station 08LF051,
# Thompson River, as exported from the WSC portal. The raw export carries one
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

thompson <- data.frame(
  date = as.Date(substr(raw$Date..UTC.[keep], 1, 10)),
  h = as.numeric(raw$Mean.Gauge.Height[keep]),
  q = as.numeric(raw$Discharge[keep]),
  uq = as.numeric(raw$Uncertainty[keep])
)
thompson <- thompson[order(thompson$h), ]
rownames(thompson) <- NULL

stopifnot(
  nrow(thompson) > 50,
  all(thompson$q > 0),
  !any(is.na(thompson$h)),
  !any(is.na(thompson$q))
)

save(thompson, file = file.path("data", "thompson.rda"), compress = "bzip2")
cat("wrote data/thompson.rda:", nrow(thompson), "gaugings\n")
