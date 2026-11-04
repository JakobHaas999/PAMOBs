### Analysis: PEMtree covariate study (split on z) ----------------------
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "pemtree_covariate_main_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

results <- load_results(run_id)
mark_valid(results)

### Checks
# Report that
sum(!is.na(results$error))
invalid_jobs <- results[valid == FALSE, .(job.id, scenario, split_type, n, delta, censoring_rate, n_events)]
invalid_jobs

# Observed censoring
censoring_check <- check_censoring(results, by = c("scenario", "split_type", "delta", "censoring_rate"))
censoring_check

### Type-I-Error rate
summary_0 <- summarise_rejection(
  results[delta == 0 & valid == TRUE],
  by = c("scenario", "split_type", "n", "censoring_rate")
)
print(summary_0)

p_type1 <- plot_rejection_by_censoring(
  summary_0,
  group = "n",
  facet = facet_thesis(rows = "split_type", cols = "scenario")
)
print(p_type1)

# Plot p value ecdf
censoring_plot <- 0.5
p_ecdf <- plot_pvalue_ecdf(
  results[delta == 0 & valid == TRUE & censoring_rate == censoring_plot],
  group = "n",
  facet = facet_thesis(rows = "split_type", cols = "scenario"),
  title = "Distribution of root-test p-values under H0",
  subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "%")
)
print(p_ecdf)

### Power analysis
summary_power <- summarise_rejection(
  results[delta > 0 & valid == TRUE],
  by = c("scenario", "split_type", "n", "censoring_rate", "delta")
)
print(summary_power)

# One plot per delta
deltas <- c(0.25, 0.5, 0.75, 1)
power_plots <- lapply(deltas, function(d) {
  plot_rejection_by_censoring(
    summary_power[abs(delta - d) < 1e-8],
    group = "n",
    facet = facet_thesis(rows = "split_type", cols = "scenario"),
    y_max = 1,
    ref_line = NULL,
    title = "Power to detect effect modification of x by z",
    subtitle = sprintf(
      "delta = %.2f | Hazard ratio for x: %.2f (z <= 0) vs. %.2f (z > 0)",
      d, exp(0.6), exp(0.6 + d)
    )
  )
})
names(power_plots) <- paste0("power_delta_", deltas)
power_plots$power_delta_0.5

# Power curves over delta (delta = 0 is the type I error)
power_curve_data <- rbind(copy(summary_0)[, delta := 0], summary_power, use.names = TRUE)
power_curves <- lapply(c("cubic", "weibull"), function(s) {
  plot_power_curve(
    power_curve_data[scenario == s],
    facet = facet_thesis(rows = "split_type", cols = "censoring_rate"),
    x_breaks = seq(0, 1, by = 0.25),
    title = paste("Power to detect effect modification of x by z |", labels_scenario[[s]]),
    subtitle = "Effect of x: 0.6 for z <= 0 and 0.6 + delta for z > 0"
  )
})
names(power_curves) <- paste0("power_curves_", c("cubic", "weibull"))
power_curves$power_curves_weibull

### Split point (numeric z, true split at 0)
tau_hat_data <- results[split_type == "numeric" & valid & split == TRUE]
summary_split_point <- summarise_split_points(
  tau_hat_data,
  truth = 0, within = 0.1,
  by = c("scenario", "n", "censoring_rate", "delta")
)
summary_split_point

# Plot split points as boxplot
censoring_plot <- 0.5
p_box <- plot_split_points_box(
  tau_hat_data[censoring_rate == censoring_plot],
  truth = 0,
  facet = facet_thesis(cols = "scenario"),
  y_limits = c(-1, 1), y_breaks = seq(-1, 1, by = 0.5),
  title = "Estimated split points (numeric z, true split at 0)",
  subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "% | only replications with a split")
)
p_box

# Histogram
n_plot <- 1000
p_hist <- plot_split_points_hist(
  tau_hat_data[n == n_plot & censoring_rate == censoring_plot],
  truth = 0,
  breaks = seq(-1, 1, by = 0.05),
  facet = facet_thesis(rows = "delta", cols = "scenario", scales = "free_y"),
  x_limits = c(-1, 1),
  title = "Distribution of estimated split points (numeric z)",
  subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%")
)
p_hist

### Save tables and plots ----------------------------------------------
save_tables(
  list(
    invalid_jobs = invalid_jobs,
    censoring_check = censoring_check,
    summary_h0 = summary_0,
    summary_power = summary_power,
    summary_split_point = summary_split_point
  ),
  results_dir
)

save_plot(p_type1, "type1_error", plot_dir)
save_plot(p_ecdf, paste0("p_value_ecdf_cens", 100 * censoring_plot), plot_dir)
for (nm in names(power_plots)) save_plot(power_plots[[nm]], nm, plot_dir)
for (nm in names(power_curves)) save_plot(power_curves[[nm]], nm, plot_dir, width = 10, height = 6.5)
save_plot(p_box, paste0("split_points_boxplot_cens", 100 * censoring_plot), plot_dir)
save_plot(p_hist, paste0("split_points_hist_n", n_plot, "_cens", 100 * censoring_plot), plot_dir, height = 9)
