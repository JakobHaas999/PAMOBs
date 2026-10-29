### setup
source("setup.R")
registry_dir <- "simulation_batchtools/registries/pilot/time_changing_effect_pilot_v1"
results_dir <- "simulation_batchtools/results/pilot/time_changing_effect_pilot_v1"

reg <- batchtools::loadRegistry(
  file.dir = registry_dir, writeable = FALSE,
  conf.file = NA, make.default = FALSE
)
job_design <- readRDS(file.path(reg$file.dir, "job_design.rds"))

done_ids <- batchtools::findDone(reg = reg)
collected <- batchtools::reduceResultsDataTable(
  ids = done_ids,
  reg = reg
)

metrics <- rbindlist(collected$result, fill = TRUE)
metrics[, job.id := collected$job.id]
metrics <- merge(as.data.table(job_design), metrics, by = "job.id", all.x = TRUE)

metrics[, completed := job.id %in% done_ids$job.id]

issues <- metrics[
  !completed | !(success %in% TRUE),
  .(job.id, scenario, seed, completed, error)
]

if (nrow(issues)) {
  print(issues)
  warning("Not all jobs completed successfully.", call. = FALSE)
}

# Create a summary table
scenario_cols <- c("scenario", "delta", "censoring_rate", "n", "df")
# 1) Null scenario
summary_table_0 <- metrics[delta == 0,
  {
    ok <- success %in% TRUE
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    power <- if (n_ok) n_splits / n_ok else NA_real_
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    list(
      n_planned = .N,
      n_success = n_ok,
      n_splits = n_splits,
      power = power,
      ci_lower = ci[1],
      ci_upper = ci[2]
    )
  },
  by = scenario_cols
]
print(summary_table_0)

# 2) Level shift scenario
summary_table_shift <- metrics[delta != 0,
  {
    ok <- success %in% TRUE
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    power <- if (n_ok) n_splits / n_ok else NA_real_
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    mean_tau_hat <- if (n_splits > 0) mean(tau_hat, na.rm = TRUE) else NA_real_
    mean_abs_tau_error <- if (n_splits) mean(abs_tau_error, na.rm = TRUE) else NA_real_

    list(
      n_planned = .N,
      n_success = n_ok,
      n_splits = n_splits,
      power = power,
      ci_lower = ci[1],
      ci_upper = ci[2],
      mean_tau_hat = mean_tau_hat,
      mean_abs_tau_error = mean_abs_tau_error
    )
  },
  by = scenario_cols
]
print(summary_table_shift)
