### Analysis: PAM tree (REML) vs. PEM tree vs. Weibull AFT tree for feature interaction (split on z)
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "pamtree_covariate_scenarios_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

res <- load_results(run_id)
mark_valid(res)

## Checks
invalid_jobs <- res[, .(n_invalid = sum(valid == FALSE), share_invalid = mean(valid == FALSE)),
  by = c("scenario", "split_type", "n", "delta", "censoring_rate", "model")
][n_invalid > 0]
invalid_jobs
censoring_check <- check_censoring(res, by = c("scenario", "split_type", "n", "delta", "censoring_rate", "model"))
censoring_check

## Type-I error (delta = 0)
summary_h0 <- summarise_rejection(
  res[delta == 0 & valid == TRUE],
  by = c("scenario", "split_type", "delta", "censoring_rate", "n", "model"), split_col = "split"
)
summary_h0

p_type1_binary <- plot_rejection_by_censoring(
  summary_h0[split_type == "binary"],
  group = "n",
  facet = facet_thesis(rows = "model", cols = "scenario"),
  title = "Binary z"
)
p_type1_binary

p_type1_numeric <- plot_rejection_by_censoring(
  summary_h0[split_type == "numeric"],
  group = "n",
  facet = facet_thesis(rows = "model", cols = "scenario"),
  title = "Numeric z"
)
p_type1_numeric

# Plot p values
censoring_plot <- 0.5
p_ecdf_binary <- plot_pvalue_ecdf(
  data = res[delta == 0 & valid == TRUE & censoring_rate == censoring_plot & split_type == "binary"],
  group = "n",
  facet = facet_thesis(rows = "model", cols = "scenario"),
  title = "Distribution of root p-values under H0",
  subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "%", " | binary split variable z")
)
p_ecdf_binary
p_ecdf_numeric <- plot_pvalue_ecdf(
  data = res[delta == 0 & valid == TRUE & censoring_rate == censoring_plot & split_type == "numeric"],
  group = "n",
  facet = facet_thesis(rows = "model", cols = "scenario"),
  title = "Distribution of root p-values under H0",
  subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "%", " | numeric split variable z")
)
p_ecdf_numeric

## 2. Power analysis
summary_power <- summarise_rejection(
  res[delta > 0 & valid == TRUE],
  by = c("scenario", "split_type", "n", "censoring_rate", "delta", "model")
)
summary_power

# Power curves over delta (delta = 0 is the type I error)
power_curve_data <- rbind(copy(summary_h0)[, delta := 0], summary_power, use.names = TRUE)
power_curves_binary <- lapply(c("cubic", "weibull"), function(s) {
  data <- power_curve_data[scenario == s & split_type == "binary"]
  plot_power_curve(
    data,
    facet = facet_thesis(rows = "model", cols = "censoring_rate"),
    x_breaks = seq(0, 1, by = 0.25),
    title = paste("Power to detect effect modification of x by z |", labels_scenario[[s]], "| Binary split variable z")
  )
})
power_curves_numeric <- lapply(c("cubic", "weibull"), function(s) {
  data <- power_curve_data[scenario == s & split_type == "numeric"]
  plot_power_curve(
    data,
    facet = facet_thesis(rows = "model", cols = "censoring_rate"),
    x_breaks = seq(0, 1, by = 0.25),
    title = paste("Power to detect effect modification of x by z |", labels_scenario[[s]], "| Numeric split variable z")
  )
})
power_curves_names <- paste0("power_curves_", c("cubic", "weibull"))
names(power_curves_binary) <- power_curves_names
names(power_curves_numeric) <- power_curves_names

### Split point (numeric z, true split at 0)
tau_hat_data <- res[split_type == "numeric" & valid == TRUE & split == TRUE]
summary_split_point <- summarise_split_points(
  tau_hat_data,
  truth = 0, within = 0.1, by = c("scenario", "n", "censoring_rate", "delta", "model")
)
summary_split_point

# Plot split points as boxplot
censoring_plot <- 0.5
p_box <- plot_split_points_box(
  tau_hat_data[delta > 0 & censoring_rate == censoring_plot],
  truth = 0,
  facet = facet_thesis(rows = "model", cols = "scenario"),
  y_limits = c(-1, 1), y_breaks = seq(-1, 1, by = 0.5),
  title = "Estimated split points (numeric z, true split at 0)",
  subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "% | only replications with a split")
)
p_box

# Histogram
n_plot <- 500 # largest n in this run
p_hist_pem <- plot_split_points_hist(
  tau_hat_data[delta > 0 & n == n_plot & censoring_rate == censoring_plot & model == "pemtree"],
  truth = 0,
  breaks = seq(-1, 1, by = 0.05),
  facet = facet_thesis(rows = "delta", cols = "scenario", scales = "free_y"),
  x_limits = c(-1, 1),
  title = "Distribution of estimated split points (numeric z)",
  subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%", "| Model = PEM tree")
)
p_hist_pem
p_hist_pam <- plot_split_points_hist(
  tau_hat_data[delta > 0 & n == n_plot & censoring_rate == censoring_plot & model == "pamtree"],
  truth = 0,
  breaks = seq(-1, 1, by = 0.05),
  facet = facet_thesis(rows = "delta", cols = "scenario", scales = "free_y"),
  x_limits = c(-1, 1),
  title = "Distribution of estimated split points (numeric z)",
  subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%", "| Model = PAM tree (REML)")
)
p_hist_pam
p_hist_aft <- plot_split_points_hist(
  tau_hat_data[delta > 0 & n == n_plot & censoring_rate == censoring_plot & model == "weibull_aft"],
  truth = 0,
  breaks = seq(-1, 1, by = 0.05),
  facet = facet_thesis(rows = "delta", cols = "scenario", scales = "free_y"),
  x_limits = c(-1, 1),
  title = "Distribution of estimated split points (numeric z)",
  subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%", "| Model = Weibull AFT tree")
)
p_hist_aft

### Save tables and plots ----------------------------------------------
save_tables(
  list(
    invalid_jobs = invalid_jobs,
    censoring_check = censoring_check,
    summary_h0 = summary_h0,
    summary_power = summary_power,
    summary_split_point = summary_split_point
  ),
  results_dir
)

save_plot(p_type1_binary, "type1_error_binary", plot_dir)
save_plot(p_type1_numeric, "type1_error_numeric", plot_dir)
save_plot(p_ecdf_binary, paste0("p_value_ecdf_binary_cens", 100 * censoring_plot), plot_dir)
save_plot(p_ecdf_numeric, paste0("p_value_ecdf_numeric_cens", 100 * censoring_plot), plot_dir)
for (nm in names(power_curves_binary)) {
  save_plot(power_curves_binary[[nm]], paste0(nm, "_binary"), plot_dir, width = 10, height = 6.5)
}
for (nm in names(power_curves_numeric)) {
  save_plot(power_curves_numeric[[nm]], paste0(nm, "_numeric"), plot_dir, width = 10, height = 6.5)
}
save_plot(p_box, paste0("split_points_boxplot_cens", 100 * censoring_plot), plot_dir)
save_plot(p_hist_pem, paste0("split_points_hist_pemtree_n", n_plot, "_cens", 100 * censoring_plot), plot_dir, height = 9)
save_plot(p_hist_aft, paste0("split_points_hist_weibull_aft_n", n_plot, "_cens", 100 * censoring_plot), plot_dir, height = 9)
