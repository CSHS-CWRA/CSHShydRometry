# Cross-check: CSHShydRometry against bdrc on the Thompson gaugings.
#
# bdrc (Hrafnkelsson et al., CRAN) fits the same power law, Q = a (h - c)^b,
# by a Bayesian hierarchical model. Its estimates and intervals should be
# close to ours where the models match:
#
#   bdrc::plm0()  constant scatter            ~ rc_powerlaw(wts = "none")
#   bdrc::plm()   scatter varying with stage  ~ rc_powerlaw(wts = "power"),
#                 roughly: ours varies with the flow, not the stage
#
# bdrc is not a dependency of CSHShydRometry; this folder is excluded from
# the package build. Run from the package root, with bdrc installed:
#
#   Rscript crosscheck/bdrc/thompson.R

pkgload::load_all(quiet = TRUE)
library(bdrc)

d <- as.data.frame(thompson[, c("discharge", "stage")])
at <- c(1, 4, 8)
set.seed(1)

ours <- list(
  none = rc_powerlaw(discharge, stage, data = d),
  power = rc_powerlaw(discharge, stage, data = d, wts = "power")
)
theirs <- list(
  none = plm0(discharge ~ stage, data = d),
  power = plm(discharge ~ stage, data = d)
)

for (nm in names(ours)) {
  cat("\n==", nm, "==\n")
  cat("\nCSHShydRometry coefficients:\n")
  print(coef(ours[[nm]]))
  cat("\nbdrc parameter summary:\n")
  print(summary(theirs[[nm]]))
  cat("\nCSHShydRometry, 95% limits:\n")
  print(predict(ours[[nm]], new_stage = at, conflev = 0.95, predlev = 0.95))
  cat("\nbdrc, 95% posterior predictive intervals:\n")
  print(predict(theirs[[nm]], newdata = at))
}
