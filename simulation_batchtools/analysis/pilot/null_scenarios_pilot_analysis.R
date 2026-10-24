registry_dir <- "simulation_batchtools/registries/pilot/null_scenario_pilot_v1"
results_dir <- "simulation_batchtools/results/pilot/null_scenario_pilot_v1"

reg <- batchtools::loadRegistry(
  file.dir = registry_dir, writeable = FALSE,
  conf.file = NA, make.default = FALSE
)
job_design <- readRDS(file.path(reg$file.dir, "job_design.rds"))

metrics <- data.table::data.table(
  job.id = integer(), success = logical(), error = character(),
  split_detected = logical(), n_terminal_nodes = integer(),
  tau_hat = numeric(), abs_tau_error = numeric(),
  observed_censoring = numeric()
)
done_ids <- batchtools::findDone(reg = reg)

if (nrow(done_ids)) {
  collected <- batchtools::reduceResultsDataTable(ids = done_ids, reg = reg)
  collected_metrics <- data.table::rbindlist(collected$result, fill = TRUE)
  collected_metrics[, job.id := collected$job.id]
  metrics <- data.table::rbindlist(list(metrics, collected_metrics), fill = TRUE)
}

metrics <- merge(
  job_design, metrics,
  by = "job.id", all.x = TRUE
)
metrics <- data.table::as.data.table(metrics)
metrics[, completed := job.id %in% done_ids$job.id]

issues <- metrics[!(success %in% TRUE), .(
  job.id, covariate_setting, scenario, rep, seed, completed, error
)]
if (nrow(issues)) {
  print(issues)
  warning("Not all jobs completed successfuly", call. = FALSE)
}

scenario_columns <- setdiff(colnames(job_design), c("job.id", "rep", "seed"))
summary_table <- metrics[,
  {
    ok <- success %in% TRUE
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])

    rejection_rate <- if (n_ok) n_splits / n_ok else NA_real_

    ci <- if (n_ok) {
      binom.test(n_splits, n_ok, p = 0.05)$conf.int
    } else {
      c(NA_real_, NA_real_)
    }

    list(
      n_planned = .N,
      n_completed = sum(completed),
      n_success = sum(ok),
      n_failed = sum(completed & !ok),
      n_missing = sum(!completed),
      n_splits = n_splits,
      rejection_rate = rejection_rate,
      ci_lower = ci[1],
      ci_upper = ci[2],
      mean_observed_censoring = if (n_ok) mean(observed_censoring[ok]) else NA_real_
    )
  },
  by = scenario_columns
]
summary_table <- summary_table[, .(
  n, censoring_rate, covariate_setting, scenario,
  n_planned, n_success, n_failed, n_missing,
  rejection_rate, ci_lower, ci_upper, mean_observed_censoring
)]

print(summary_table)

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(metrics, file.path(results_dir, "job_results.rds"))
data.table::fwrite(summary_table, file.path(results_dir, "summary.csv"))
