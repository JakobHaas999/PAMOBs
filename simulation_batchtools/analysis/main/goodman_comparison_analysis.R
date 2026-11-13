### Analysis of pemtree with goodman et al (2011)
### setup
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "goodman_comparison_v1"
resuts_dir <- results_path(run_id)
plot_dir <- file.path(resuts_dir, "figures")

res <- load_results(run_id)
mark_valid(res)

## 1. Checks
invalid_jobs <- res[, .(n_invalid = sum(!valid), share_invalid = mean(!valid)),
  by = c("n", "delta", "censoring_rate", "model")
][n_invalid > 0]
invalid_jobs
censoring_check <- check_censoring(res, by = c("n", "delta", "censoring_rate", "model"))
censoring_check

## 2. Type-I error (delta = 0)
summary_h0 <- summarise_rejection(
  res[delta == 0 & valid == TRUE],
  by = c("n", "model", "delta", "censoring_rate")
)
summary_h0

p_type1 <- plot_rejection_by_censoring(
  summary_h0,
  group = "n", facet = facet_thesis(cols = "model")
)

# Plot p values of PEM
censoring_plot <- 0.5
p_ecdf <- plot_pvalue_ecdf(
  res[delta == 0 & censoring_rate == censoring_plot & valid == TRUE & model == "pemtree"],
  group = "n",
  title = "Distribution of p values under H0",
  subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "%")
)

## 3. Power analysis
summary_power <- summarise_rejection(
  res[delta > 0 & valid == TRUE],
  by = c("n", "model", "delta", "censoring_rate")
)
summary_power

# Power curves over delta (delta = 0 is the type-I error)
power_curve_data <- rbind(
  copy(summary_h0)[, delta := 0],
  summary_power,
  use.names = TRUE
)

p_power <- plot_power_curve(
  power_curve_data,
  facet = facet_thesis(rows = "model", cols = "censoring_rate"),
  x_breaks = sort(unique(power_curve_data$delta))
)

## 4- Split point (trur split at tau = 3)
tau_hat_data <- res[split == TRUE & valid == TRUE]
summary_split_point <- summarise_split_points(
  tau_hat_data,
  truth = 3,
  within = 0.5,
  by = c("n", "delta", "censoring_rate", "model")
)

# Plot split points as boxplot
p_box <- plot_split_points_box(
  tau_hat_data[delta > 0],
  truth = 3,
  facet = facet_thesis(rows = "model", cols = "censoring_rate"),
  y_limits = c(0, 10), y_breaks = seq(0, 10, by = 1),
  title = "Estimated split points (true split at 3)",
  subtitle = "Only replications with a split"
)

# Histogram
n_plot <- 500
p_hist <- plot_split_points_hist(
  tau_hat_data[delta > 0 & n == n_plot & censoring_rate == censoring_plot],
  truth = 3,
  breaks = seq(-0.05, 10.05, by = 0.25),
  facet = facet_thesis(rows = "model", cols = "delta", scales = "free_y"),
  x_limits = c(0, 10),
  title = "Distribution of estimated split points (true split at 3)",
  subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%")
)
p_hist

## 5. Save tables and plots
save_tables(
  list(
    invalid_jobs = invalid_jobs,
    censoring_check = censoring_check,
    summary_h0 = summary_h0,
    summary_power = summary_power,
    summary_split_point = summary_split_point
  ),
  resuts_dir
)

save_plot(p_type1, "type1_error", plot_dir, width = 8, height = 4.5)
save_plot(p_ecdf, paste0("p_value_ecdf_pemtree_cens", 100 * censoring_plot), plot_dir)
save_plot(p_power, "power_curves", plot_dir, width = 10, height = 6.5)
save_plot(p_box, "split_points_boxplot", plot_dir, width = 10, height = 6.5)
save_plot(p_hist, paste0("split_points_hist_n", n_plot, "_cens", 100 * censoring_plot), plot_dir, width = 10, height = 5)
