### setup
source("setup.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemtree_splines.R")
source("R/pamob.R")
# Simulation helper
source("simulation/R/run_covariate_instability.R")

# Parameters
log_hazard <- function(t) {
  -3 + 0.15 * t - 0.01 * t^2
}

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -1, max = 1),
  x3 = list(distfun = rbinom, prob = 0.5, size = 1)
)

design <- expand.grid(
  n = c(1000, 2000),
  censoring_rate = 0.3,
  rep = seq_len(50),
  admin_time = 10,
  sim_interval = 0.05,
  fit_interval = 0.1,
  KEEP.OUT.ATTRS = FALSE
)

sim_formula <- ~ log_hazard(t) + 0.7 * x1 + 0.4 * x2
fit_formula <- Surv(time, status) ~ poly(tend, degree = 2, raw = TRUE) + x1 + x2

# Run pilot simulation
results_pamtree <- run_covariate_instability(
  log_hazard = log_hazard,
  design = design,
  sim_formula = sim_formula,
  fit_formula = fit_formula,
  fit_fun = pemtree_splines,
  covariate_spec = covariate_spec,
  alpha = 0.05,
  maxdepth = 2
)
results_pemob <- run_covariate_instability(
  log_hazard = log_hazard,
  design = design,
  sim_formula = sim_formula,
  fit_formula = fit_formula,
  fit_fun = pemob,
  covariate_spec = covariate_spec,
  alpha = 0.05,
  maxdepth = 2
)
results <- rbind.data.frame(results_pamtree, results_pemob)
# Save results
saveRDS(results, "simulation/results/pilot/07_compare_pemob_pemtree_pilot.rds")

# Analysis of results
