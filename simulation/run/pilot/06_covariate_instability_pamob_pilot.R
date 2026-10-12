### setup
source("setup.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pamob.R")
# Simulation helper
source("simulation/R/run_covariate_instability.R")


log_hazard <- function(t) {
  -1.1 * sin(pi / 4) * t
}

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = 0, max = 3),
  x3 = list(distfun = rbinom, prob = 0.5, size = 1)
)

# admin_time <- 10
# cut <- seq(0, admin_time, by = 0.05)

design <- expand.grid(
  n = c(1000, 2000),
  delta = c(0, 0.4, 0.8),
  censoring_rate = 0.3,
  rep = seq_len(10),
  admin_time = 10,
  sim_interval = 0.05,
  fit_interval = 0.1,
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 8000 + seq_len(nrow(design))

sim_formula <- function(d) {
  reformulate(
    termlabels = c(
      "log_hazard(t)",
      "0.7 * x1",
      "0.4 * x2",
      sprintf("(%s) * (x3 == 1) * x1", d$delta)
    ),
    env = environment()
  )
}
fit_formula <- Surv(time, status) ~ s(tend, bs = "ps", k = 8) + x1 + x2 + offset(offset) | x3

results <- run_covariate_instability(
  log_hazard = log_hazard,
  design = design,
  sim_formula = sim_formula,
  fit_formula = fit_formula,
  fit_fun = pamob,
  covariate_spec = covariate_spec,
  min_events = 20,
  alpha = 0.05,
  maxdepth = 2,
  vcov = "opg"
)
saveRDS(results, "simulation/results/pilot/06_covar_instability_pamob_pilot.rds")

# Inspect results
data.table::setDT(results)

# Null scenario
results[success == TRUE & delta == 0, .(
  alpha_hat = mean(split_detected),
  mcse = sqrt(mean(split_detected) * (1 - mean(split_detected)) / .N)
), by = .(n)]
# Power
results[success == TRUE & delta != 0, .(
  power_hat = mean(split_detected),
  mcse = sqrt(mean(split_detected) * (1 - mean(split_detected)) / .N)
), by = .(n, delta)]
