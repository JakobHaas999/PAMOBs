source("setup.R")
## Analysis of results
results_dir <- "simulation_custom/results"
results <- readRDS(file.path(results_dir, "weibul_level_shift.rds"))

# Create a summary table
results <- as.data.table(results)
summary <- results[,
  {
    ok <- success %in% TRUE &
      is.finite(root_p_value)
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    power <- if (n_ok) n_splits / n_ok else NA_real_
    ci = if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)
    
    has_tau <- ok & (split_detected %in% TRUE) & is.finite(tau_hat)
    if (any(has_tau)) {
      tau_error <- mean(abs(tau_hat[has_tau] - tau[has_tau]))
      mean_tau_hat <- mean(tau_hat[has_tau])
    } else {
      tau_error <- NA_real_
      mean_tau_hat <- NA_real_
    }

    list(
      n_planned = .N, n_success = n_ok, n_split = n_splits,
      power = if (n_ok) n_splits / n_ok else NA_real_,
      ci_lower = ci[1], ci_upper = ci[2],
      n_tau = sum(has_tau),
      mean_tau_hat = if (any(has_tau)) mean(tau_hat[has_tau]) else NA_real_,
      tau_mae = if (any(has_tau)) mean(abs(tau_hat[has_tau] - tau[has_tau])) else NA_real_
    )
    
}, by = c("model", "hazard_ratio")]
summary