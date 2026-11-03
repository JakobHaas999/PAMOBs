source("setup.R")
### Analysis of results
results_dir <- "simulation_custom/results"
results <- readRDS(file.path(results_dir, "pamtree_time_changing_effect_main.rds"))

sim_results <- results$results
sim_results$tau_true <- 4

# Create summary table for null scenario
summary_0 <- sim_results[delta == 0, {
  ok <- success %in% TRUE & is.finite(root_p_value)
  n_ok <- sum(ok)
  n_splits <- sum(split_detected[ok])
  split_rate <- if (n_ok) n_splits / n_ok else NA_real_
  ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

  list(
    n_planned = .N, n_success = n_ok, n_split = n_splits,
    split_rate = split_rate,
    ci_lower = ci[1], ci_upper = ci[2]
  )
}]
print(summary_0)

# Create summary table for level shift scenario
summary_level_shift <- sim_results[delta != 0, {
  ok <- success %in% TRUE & is.finite(root_p_value)
  n_ok <- sum(ok)
  n_splits <- sum(split_detected[ok])
  split_rate <- if(n_ok) n_splits / n_ok else NA_real_
  ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

  has_tau <- ok & (split_detected %in% TRUE) & is.finite(tau_hat)
  if (any(has_tau)) {
    tau_error = mean(abs(tau_hat[has_tau] - tau_true[has_tau]))
    mean_tau_hat = mean(tau_hat[has_tau])
  } else {
    tau_error <- NA_real_
    mean_tau_hat <- NA_real_
  }

  list(
    n_planned = .N, n_success = n_ok, n_split = n_splits,
    power = split_rate, ci_lower = ci[1], ci_upper = ci[2],
    n_tau = sum(has_tau),
    mean_tau_hat = mean_tau_hat,
    tau_error = tau_error
  )
}, by = "delta"]
print(summary_level_shift)

