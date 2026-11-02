library("data.table")
library("batchtools")

reg <- suppressWarnings(
  loadRegistry("simulation_batchtools/registries/pilot/pemtree_covariate_pilot_v1", work.dir = getwd())
)

results <- merge(
  unwrap(getJobPars()),
  unwrap(reduceResultsDataTable()),
  by = "job.id"
)

nrow(results[!is.na(error)])
nrow(results[is.na(split)])
obs_cens <- results[, .(mean_observed_censoring = mean(1 - n_events / n)), by = censoring_rate]
obs_cens

summary_table <- results[, .(
  rate = mean(split),
  mc_se = sqrt(mean(split) * (1 - mean(split)) / .N),
  mean_tau_hat = mean(split_point, na.rm = TRUE),
  sd_tau_hat = sd(split_point, na.rm = TRUE)
), by = .(scenario, split_type, delta, censoring_rate, n)]
summary_table
