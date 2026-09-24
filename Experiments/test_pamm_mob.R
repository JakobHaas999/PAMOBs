source("setup.R")
source("R/simulation_pam.R")

######## Experiment:
### We want to find a structural break in the baseline hazard
### using a fixed spline basis, no penalization
### We want to find the break through model based partitioning

## Simulate some data
set.seed(12)
log_baseline <- function(t) {
  -2 + 0.15 * t + ifelse(t > 4, 3, 0)
}
formula <- ~ log_baseline(t)
cut <- seq(0, 8, by = 0.25)
surv_data <- sim_pam(
  n = 10000,
  formula = formula,
  covariate_spec = NULL,
  baseline = log_baseline,
  cut = cut
)

ped <- as_ped(Surv(time, status) ~ 1, data = surv_data, cut = cut)
