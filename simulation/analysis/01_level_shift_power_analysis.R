##### Analysis of 01 level shift power

source("setup.R")

## Load results data

# Bind blocks files to one file
run_id <- "level_shift_power_robust_v1"
results_folder <- file.path("simulation/results/main", run_id, "blocks")
blocks <- list.files(results_folder, pattern = "^block_[0-9]+\\.rds$")
results <- rbindlist(lapply(blocks, function(b) {
  as.data.table(readRDS(file.path(results_folder, b)))
}))

# Check results
assert_true(nrow(results) == 12000)
assert_true(all(results$success))
assert_true(!anyDuplicated(results$design_id))

# Summary of results
results_summary <- results[
  success == TRUE,
  .(
    n_sim = .N,
    power = mean(split_detected),
    cond_mean_tau_hat = mean(tau_hat, na.rm = TRUE),
    cond_mean_tau_error = mean(tau_error, na.rm = TRUE),
    cond_bias_tau = mean(tau_hat - tau, na.rm = TRUE),
    cond_rmse_tau = sqrt(mean((tau_hat - tau)^2, na.rm = TRUE))
  ),
  by = c("delta", "n", "censoring_rate")
][, mcse_power := sqrt(power * (1 - power) / n_sim)][]


# Power summary table
power_summary <- results_summary[
  ,
  .(
    delta, n, censoring_rate, power,
    power_lower = pmax(0, power - 1.96 * mcse_power),
    power_upper = pmin(1, power + 1.96 * mcse_power)
  )
]

# Power curves
ggplot(power_summary, aes(x = delta, y = power, colour = factor(n), group = factor(n))) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey80") +
  geom_ribbon(aes(ymin = power_lower, ymax = power_upper, fill = factor(n)),
    alpha = 0.15, colour = NA
  ) +
  geom_line() +
  geom_point() +
  facet_wrap(~censoring_rate) +
  scale_x_continuous(breaks = c(0, unique(power_summary$delta)), limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(x = expression(delta), y = "Empirical Power", color = "Sample size", fill = "Sample size") +
  theme_bw()

# Precision of split point
ggplot(
  results[split_detected == TRUE],
  aes(x = factor(delta), y = tau_hat, fill = factor(n))
) +
  geom_hline(yintercept = 4, linetype = "dashed") +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~censoring_rate) +
  labs(x = expression(delta), y = expression(hat(tau)), fill = "Sample size") +
  theme_bw()
