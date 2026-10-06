### source relevant files
## packages
source("setup.R")

## functions
source("R/pemtree.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")

## simulation helpers
source("simulation/R/fit_models.R")
source("simulation/R/run_one.R")

## simulation scenario
source("simulation/scenarios/01_level_shift.R")


### Pilot run
# Create design
design <- expand.grid(
  n = 500,
  delta = c(0, 0.5, 1),
  tau = 4,
  beta0 = -2,
  censoring_rate = 0.3,
  rep = seq_len(10),
  ped_interval = 0.1,
  alpha = 0.05,
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 1000 + seq_len(nrow(design))

results <- do.call(
  rbind,
  lapply(seq_len(nrow(design)), function(i) {
    do.call(
      run_one_level_shift,
      as.list(design[i, ])
    )
  })
)

results_summary <- results |>
  group_by(delta) |>
  summarise(
    n_runs = n(),
    n_errors = sum(!success),
    split_rate = if (any(success)) {
      mean(split_detected[success])
    } else {
      NA_real_
    },
    mean_tau_hat = if (all(is.na(tau_hat))) {
      NA_real_
    } else {
      mean(tau_hat, na.rm = TRUE)
    },
    mean_tau_error = if (all(is.na(tau_error))) {
      NA_real_
    } else {
      mean(tau_error, na.rm = TRUE)
    }
  )

cat("Summary of results:\n\n")
print(results_summary)

## Save results
saveRDS(
  object = results,
  file = "simulation/results/pilot/01_level_shift_raw.rds"
)
