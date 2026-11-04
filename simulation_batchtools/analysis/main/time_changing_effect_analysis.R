### Analysis: PEMtree time-varying effect of x (split on time) ---------
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "time_changing_effect_main_v2"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

results <- load_results(run_id)
mark_valid(results)

### Checks
sum(!is.na(results$error))
invalid_jobs <- results[valid == FALSE, .(job.id, scenario, n, delta, censoring_rate, n_events)]
invalid_jobs

# Observed censoring
censoring_check <- check_censoring(results, by = c("scenario", "delta", "n", "censoring_rate"))
censoring_check

### Type-I-Error rate
summary_0 <- summarise_rejection(
  results[delta == 0 & valid == TRUE],
  by = c("scenario", "n", "censoring_rate")
)
summary_0

p_type1 <- plot_rejection_by_censoring(
  summary_0,
  group = "n",
  facet = facet_thesis_wrap("scenario")
)
p_type1

### Power analysis
summary_power <- summarise_rejection(
  results[delta > 0 & valid == TRUE],
  by = c("scenario", "n", "censoring_rate", "delta")
)
summary_power

# One plot per delta
deltas <- c(0.25, 0.5, 0.75, 1)
power_plots <- lapply(deltas, function(d) {
  plot_rejection_by_censoring(
    summary_power[abs(delta - d) < 1e-8],
    group = "n",
    facet = facet_thesis_wrap("scenario"),
    y_max = 1,
    ref_line = NULL,
    title = "Power to detect a time-varying effect of x",
    subtitle = sprintf(
      "delta = %.2f | Hazard ratio for x: %.2f (t < 3) vs. %.2f (t >= 3)",
      d, exp(0.6), exp(0.6 + d)
    )
  )
})
names(power_plots) <- paste0("power_delta_", deltas)
power_plots$power_delta_0.5

### Power curves over delta (delta = 0 is the type I error)
power_curve_data <- rbind(copy(summary_0)[, delta := 0], summary_power, use.names = TRUE)
p_power_curve <- plot_power_curve(
  power_curve_data,
  facet = facet_thesis(rows = "scenario", cols = "censoring_rate"),
  x_breaks = seq(0, 1, by = 0.25),
  title = "Power to detect a time-varying effect of x",
  subtitle = "Effect of x: 0.6 before t = 3 and 0.6 + delta afterwards"
)
p_power_curve

### Split point analysis (true change at t = 3)
tau_hat_data <- results[valid == TRUE & split == TRUE & delta > 0]
summary_split_point <- summarise_split_points(
  tau_hat_data,
  truth = 3,
  by = c("scenario", "n", "censoring_rate", "delta")
)
summary_split_point

# Plot split points as boxplot
censoring_plot <- 0.5
p_box <- plot_split_points_box(
  tau_hat_data[censoring_rate == censoring_plot],
  truth = 3,
  facet = facet_thesis_wrap("scenario"),
  y_limits = c(0, 10), y_breaks = seq(0, 10, by = 2),
  title = "Estimated split points in time (true change at t = 3)",
  subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "% | only replications with a split")
)
p_box

# Histogram of split points (Weibull, n = 1000)
n_plot <- 1000
p_hist <- plot_split_points_hist(
  tau_hat_data[scenario == "weibull" & n == n_plot & censoring_rate == censoring_plot],
  truth = 3,
  breaks = seq(-0.05, 10.05, by = 0.1),
  facet = facet_thesis_wrap("delta", ncol = 1, scales = "free_y"),
  x_limits = c(0, 10), x_breaks = 0:10,
  title = "Distribution of estimated split points (Weibull baseline hazard)",
  subtitle = paste0(
    "n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "% | true change at t = 3"
  )
)
p_hist

### Coefficients of x before and after the split
coef_data <- results[valid & split == TRUE & delta > 0]
summary_coef <- coef_data[, .(
  n_splits = .N,
  mean_before = mean(x_before),
  bias_before = mean(x_before - beta),
  mean_after = mean(x_after),
  true_after = beta[1] + delta[1],
  bias_after = mean(x_after - (beta + delta)),
  rmse_before = sqrt(mean((x_before - beta)^2)),
  rmse_after = sqrt(mean((x_after - (beta + delta))^2))
), by = c("scenario", "n", "censoring_rate", "delta")][
  order(scenario, n, censoring_rate, delta)
]
summary_coef

# Boxplot of estimates vs. true values (n = 1000, 50 % censoring)
coef_long <- melt(
  coef_data[n == 1000 & censoring_rate == 0.5],
  id.vars = c("scenario", "delta", "beta"),
  measure.vars = c("x_before", "x_after"),
  variable.name = "period", value.name = "estimate"
)
coef_long[, truth := fifelse(period == "x_before", beta, beta + delta)]
coef_long[, period := factor(period,
  levels = c("x_before", "x_after"),
  labels = c("Before t = 3", "After t = 3")
)]
truth_data <- unique(coef_long[, .(scenario, delta, period, truth)])

p_coef <- ggplot(coef_long, aes(x = factor(delta), y = estimate, fill = period)) +
  geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(
    data = truth_data, aes(x = factor(delta), y = truth, group = period),
    position = position_dodge(width = 0.8), shape = 4, size = 3, stroke = 1.2,
    colour = "black", inherit.aes = FALSE
  ) +
  facet_thesis_wrap("scenario") +
  scale_fill_manual(values = c("Before t = 3" = "grey80", "After t = 3" = "#56B4E9")) +
  labs(
    title = "Estimated effect of x in the two terminal nodes",
    subtitle = "n = 1000 | Target censoring rate = 50% | only replications with a split | x = true value",
    x = expression(delta), y = "Estimated coefficient of x", fill = NULL
  ) +
  theme_thesis(discrete_x = TRUE)
p_coef

### Save tables and plots ----------------------------------------------
save_tables(
  list(
    invalid_jobs = invalid_jobs,
    censoring_check = censoring_check,
    summary_h0 = summary_0,
    summary_power = summary_power,
    summary_split_point = summary_split_point,
    summary_coef = summary_coef
  ),
  results_dir
)

save_plot(p_type1, "type1_error", plot_dir)
for (nm in names(power_plots)) save_plot(power_plots[[nm]], nm, plot_dir)
save_plot(p_power_curve, "power_curves", plot_dir, width = 10, height = 6.5)
save_plot(p_box, paste0("split_points_boxplot_cens", 100 * censoring_plot), plot_dir)
save_plot(p_hist, paste0("split_points_hist_weibull_n", n_plot, "_cens", 100 * censoring_plot), plot_dir, height = 8)
save_plot(p_coef, "coefficients_x_n1000_cens50", plot_dir)
