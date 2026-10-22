# Thompson River gaugings

Stage and discharge measurements from Water Survey of Canada station
08LF051, Thompson River near Spences Bridge, British Columbia, used in
the examples: every discharge measurement made at the station from 1994
to 2024, in date order.

## Usage

``` r
thompson
```

## Format

A tibble with 93 rows and 5 columns:

- date:

  Date of the measurement.

- stage:

  Stage, in meters (the mean gauge height for the visit).

- discharge:

  Discharge, in cubic meters per second.

- uncertainty_pct:

  Reported uncertainty of the discharge, as a percentage of it at two
  standard deviations (the "IVE method, 2-sigma value"), or `NA` where
  none was given.

- rating_table:

  The Water Survey's rating table in force at the gauging, such as
  `"11"`, or `NA` where none was recorded (all gaugings before 2020).

## Source

Environment and Climate Change Canada (2026). Historical hydrometric
data: Station 08LF051, Thompson River near Spences Bridge, BC (open
water field visits, 1994-2024) \[Data set\]. Water Survey of Canada.
<https://wateroffice.ec.gc.ca/>

Included with the permission of the Water Survey of Canada.

## Details

Reported uncertainties are available for only 19 of the gaugings, so
`uncertainty_pct` is mostly `NA`. It is a percentage of the discharge at
two standard deviations, so the standard uncertainty of a discharge, in
cubic meters per second, is `uncertainty_pct / 100 * discharge / 2`.
Fits using
[`var_spec()`](https://cshs-cwra.github.io/CSHShydRometry/reference/variance.md)
need a variance for every gauging, and so need either the subset with an
uncertainty or variances from elsewhere.

The rating at the station has shifted over the decades: gaugings from
different periods depart from a single fitted curve by several percent
in opposite directions, and the Water Survey's rating table changed
three times between 2020 and 2024. A curve fitted to all 93 gaugings is
fine for illustration, but a rating for a given period should use
gaugings from a period of stable channel conditions; see
[`vignette("fitting")`](https://cshs-cwra.github.io/CSHShydRometry/articles/fitting.md).
