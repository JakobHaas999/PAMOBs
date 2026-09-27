### source relevant files
## packages
source("setup.R")

## functions
source("R/pemtree.R")
source("R/sim_pexp_own.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")

## simulation helpers
source("simulation/R/fit_models.R")
source("simulation/R/run_one.R")

## simulation scenario
source("simulation/scenarios/01_level_shift.R")

### Pilot run
design <- expand.grid(
  n = 500,
  delta = 0,
  tau = 4,
  beta0 = -2,
  censoring_rate = 0.3,
  rep = seq_len(200),
  ped_interval = 0.1,
  alpha = 0.05,
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 10000 + seq_len(nrow(design))

results <- do.call(
  rbind,
  lapply(seq_len(nrow(design)), function(i) {
    do.call(
      run_one_level_shift,
      as.list(design[i, ])
    )
  })
)

## Summary statistics
split_rate <- mean(results$split_detected[results$success])
mcse = sqrt(split_rate * (1 - split_rate) / sum(results$success))
c(
  n_errors = sum(!results$success),
  split_rate = split_rate,
  mcse = mcse
)
results |>
  filter(success, split_detected) |>
  select(success, tau_hat)

## Save results
saveRDS(
  object = results,
  "simulation/results/pilot/02_level_shift_null_raw.rds"
)