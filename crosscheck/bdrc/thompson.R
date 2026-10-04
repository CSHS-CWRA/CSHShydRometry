# Cross-check: CSHShydRometry against bdrc on the Thompson gaugings.
#
# bdrc (Hrafnkelsson et al., CRAN) fits the same power law, Q = a (h - c)^b,
# by a Bayesian hierarchical model on the log scale. So its "constant
# scatter" model is constant on the log scale, which is scatter proportional
# to the flow. The comparable fits are:
#
#   bdrc::plm0()  constant scatter on the log scale  ~ rc_powerlaw_log(),
#                                                      rc_powerlaw(variance = "prop")
#   bdrc::plm()   log-scale scatter varying with     ~ rc_powerlaw(variance = "power"),
#                 stage                                roughly: ours varies with
#                                                      the flow, not the stage
#
# bdrc's priors are weak, so its intervals should be close to ours.
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
  powerlaw_log = rc_powerlaw_log(discharge, stage, data = d),
  powerlaw_prop = rc_powerlaw(discharge, stage, data = d, variance = "prop"),
  powerlaw_power = rc_powerlaw(discharge, stage, data = d, variance = "power")
)
theirs <- list(
  plm0 = plm0(discharge ~ stage, data = d, verbose = FALSE),
  plm = plm(discharge ~ stage, data = d, verbose = FALSE)
)

# parameters: our estimates, and bdrc's posterior medians and 95% intervals
pars <- rbind(
  t(sapply(ours, coef)),
  t(sapply(theirs, function(f) f$param_summary[c("a", "b", "c"), "median"]))
)
cat("\nParameters (bdrc: posterior medians)\n")
print(round(pars, 3))
for (nm in names(theirs)) {
  cat("\nbdrc", nm, "95% intervals\n")
  print(round(theirs[[nm]]$param_summary[c("a", "b", "c"), c("lower", "upper")], 3))
}

# limits on the curve: our confidence limits, and bdrc's credible limits for
# a (h - c)^b from its posterior draws
curve_limits <- function(f) {
  draws <- sapply(at, function(h) {
    f$a_posterior * pmax(h - f$c_posterior, 0)^f$b_posterior
  })
  t(apply(draws, 2, stats::quantile, c(0.025, 0.975)))
}
limits <- list()
for (nm in names(ours)) {
  p <- predict(ours[[nm]], new_stage = at, conflev = 0.95, predlev = 0.95)
  limits[[nm]] <- data.frame(model = nm, stage = at, fit = p$fit,
    ci_lwr = p$ci_lwr, ci_upr = p$ci_upr, pi_lwr = p$pi_lwr, pi_upr = p$pi_upr)
}
for (nm in names(theirs)) {
  f <- theirs[[nm]]
  ci <- curve_limits(f)
  pi <- predict(f, newdata = at)
  limits[[nm]] <- data.frame(model = nm, stage = at, fit = pi$median,
    ci_lwr = ci[, 1], ci_upr = ci[, 2], pi_lwr = pi$lower, pi_upr = pi$upper)
}
limits <- do.call(rbind, limits)
limits <- limits[order(limits$stage), ]
rownames(limits) <- NULL
limits[-(1:2)] <- round(limits[-(1:2)])
cat("\n95% limits (bdrc: credible limits for the curve, posterior predictive",
  "limits for a new gauging)\n")
print(limits)
