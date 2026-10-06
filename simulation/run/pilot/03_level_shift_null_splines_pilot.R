#### setup
# packages
source("setup.R")
# functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemtree.R")
source("R/pemtree_splines.R")
# simulation helpers
source("simulation/R/pemtree_helpers.R")
source("simulation/R/run_one_level_shift_splines.R")
source("simulation/scenarios/02_level_shift_splines.R")

###### Pilot run
design <- expand.grid(
  n = 1000,
  delta = 0,
  tau = 5,
  covariate_setting = c("none", "included"),
  df = c(4, 5, 6),
  censoring_rate = 0.2,
  rep = seq_len(20),
  admin_time = 10,
  sim_interval = 0.05,
  ped_interval = 0.1,
  min_events = 25,
  alpha = 0.05,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design$seed <- 300000 + seq_len(nrow(design))

basehaz_fun <- function(t) {
  -2.5 + 0.25 * t - 0.025 * t^2
}

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -1, max = 1),
  x3 = list(distfun = rbinom, size = 1, prob = 0.5)
)

covariate_effects <- c(x1 = 0.3, x2 = 0.2, x3 = 0.5)

results <- do.call(
  rbind,
  lapply(seq_len(nrow(design)), function(i) {
    args <- as.list(design[i, , drop = FALSE])
    args$basehaz_fun <- basehaz_fun

    if (args$covariate_setting == "included") {
      args$covariate_spec <- covariate_spec
      args$covariate_effects <- covariate_effects
    }

    do.call(run_one_level_shift_splines, args)
  })
)

saveRDS(
  object = results,
  "simulation/results/pilot/03_level_shift_null_splines_raw.rds"
)

## Summary statistics
results <- as.data.table(results)
results[success == TRUE, .(
  alpha_hat = mean(split_detected),
  mcse = sqrt(mean(split_detected) * (1 - mean(split_detected)) / .N)
), by = .(scenario, df)]

