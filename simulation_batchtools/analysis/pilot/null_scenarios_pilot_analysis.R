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
  job_design, metrics, by = "job.id", all.x = TRUE
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
