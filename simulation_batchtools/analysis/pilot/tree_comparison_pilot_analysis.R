registry_dir <- "simulation_batchtools/registries/pilot/tree_comparison_pilot_v1"
results_dir <- "simulation_batchtools/results/pilot/tree_comparison_pilot_v1"

reg <- batchtools::loadRegistry(
  file.dir = registry_dir,
  writeable = FALSE,
  conf.file = NA,
  make.default = FALSE
)

job_design <- readRDS(file.path(reg$file.dir, "job_design.rds"))

collected <- batchtools::reduceResultsDataTable(
  ids = batchtools::findDone(reg = reg),
  reg = reg
)

metrics <- data.table::rbindlist(collected$result, fill = TRUE)
metrics[, job.id := collected$job.id]
metrics <- merge(data.table::as.data.table(job_design), metrics, by = "job.id")

summary_table <- metrics[, {
  ok <- success %in% TRUE &
    test_available %in% TRUE &
    is.finite(root_p_value) &
    !is.na(split_detected)

  n_ok <- sum(ok)
  n_splits <- sum(split_detected[ok])

  ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

  list(
    n_completed = .N,
    n_valid = n_ok,
    n_failed = sum(success %in% FALSE),
    n_missing_test = sum(
      success %in% TRUE & !(test_available %in% TRUE)
    ),
    n_splits = n_splits,
    split_rate = if (n_ok) n_splits / n_ok else NA_real_,
    ci_lower = ci[[1]],
    ci_upper = ci[[2]]
  )
}, by = .(scenario, censoring_rate, delta, n, df, method)]

data.table::setorder(summary_table, scenario, censoring_rate, delta, method)
print(summary_table[delta == 0])
print(summary_table[delta == 0.4])
print(summary_table[delta == 0.8])

comparison <- data.table::dcast(
  summary_table,
  scenario + censoring_rate + delta + n + df ~ method,
  value.var = "split_rate"
)

print(comparison)