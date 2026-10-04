# Cross-checks

Comparisons of CSHShydRometry's fits with other software, run by hand.
They are not tests, and the other software is not a dependency: this
folder is excluded from the package build.

- `bdrc/thompson.R`: power-law fits to the Thompson gaugings, against the
  Bayesian fits of the bdrc package. Install bdrc first.

## Results

`bdrc/thompson.R`, run 2026-10-04 with bdrc 2.0.1. bdrc's priors are weak,
and its results agree closely with ours. Medians match our estimates to
within 1% (a ≈ 52.5, b = 1.88, c ≈ −1.30), and the 95% limits match to
within a few m³/s:

| Stage | Ours, `variance = "power"`: curve | bdrc `plm0()`: curve | Ours: new gauging | bdrc: new gauging |
|---|---|---|---|---|
| 1 m | 249–255 | 248–255 | 232–272 | 232–273 |
| 4 m | 1195–1222 | 1195–1221 | 1111–1306 | 1112–1306 |
| 8 m | 3410–3542 | 3408–3538 | 3191–3761 | 3197–3772 |

This independently supports the delta-method limits for `var_power()`.
nlraa's simulated limits, used before, gave 1178–1218 for the curve at 4 m.
