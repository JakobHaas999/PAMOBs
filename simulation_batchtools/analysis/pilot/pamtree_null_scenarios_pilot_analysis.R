library("data.table")
library("batchtools")

reg <- suppressWarnings(
  loadRegistry("simulation_batchtools/registries/pilot/pamtree_null_scenarios_pilot_v1", work.dir = getwd())
)

res <- reduceResultsDataTable()
res <- res[, result[[1]], by = "job.id"]

results <- merge(unwrap(getJobPars()), res, by = "job.id")

summary_h0 <- results[,
  {
    ok <- is.na(error) & is.finite(p_root)
    n_ok <- sum(ok)

    n_splits <- sum(split[ok])
    split_rate <- if (n_ok) n_splits / n_ok else NA_real_
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    list(
      n_planned = .N,
      n_completed = n_ok,
      split_rate = split_rate,
      ci_lower = ci[1],
      ci_upper = ci[2]
    )
  },
  by = c("variant", "censoring_rate")
][order(variant, censoring_rate)]
print(summary_h0)
