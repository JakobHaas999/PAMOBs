### Analysis of exposure-scaled pemtree with goodman et al (2011)
### setup
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "goodman_es_comparison_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

res <- load_results(run_id)
mark_valid(res)

## 1. Checks
invalid_jobs <- res[, .(n_invalid = sum(!valid), share_invalid = mean(!valid)),
  by = c("n", "delta", "censoring_rate", "model")
][n_invalid > 0]
invalid_jobs
censoring_check <- check_censoring(res, by = c("n", "delta", "censoring_rate", "model"))
censoring_check

## 2. Type I error (delta = 0)
summary_h0 <- summarise_rejection(
  res[delta == 0 & valid == TRUE],
  by = c("n", "model", "delta", "censoring_rate")
)
summary_h0

p_type1 <- plot_rejection_by_censoring(
  summary_h0,
  group = "n",
  facet = facet_thesis(cols = "model")
)

## 3. Power analysis
summary_power <- summarise_rejection(
  res[delta > 0 & valid == TRUE],
  by = c("n", "model", "delta", "censoring_rate")
)
summary_power

# Power curves over delta
power_curve_data <- rbind(
  copy(summary_h0)[, delta := 0],
  summary_power,
  use.names = TRUE
)

p_power <- plot_power_curve(
  power_curve_data[censoring_rate == 0.5],
  facet = facet_thesis(rows = "model"),
  x_breaks = sort(unique(power_curve_data$delta))
)
p_power

## 4. Split point (true split at tau = 3)
tau_hat_data <- res[split == TRUE & valid == TRUE]
summary_split_point <- summarise_split_points(
  tau_hat_data,
  truth = 3,
  within = 0.5,
  by = c("n", "delta", "censoring_rate", "model")
)

# Plot split points as boxplot
p_box <- plot_split_points_box(
  tau_hat_data[delta > 0 & censoring_rate == 0.5],
  truth = 3,
  facet = facet_thesis(rows = "model")
)
p_box

# Histgram
n_plot <- 500
p_hist <- plot_split_points_hist(
  tau_hat_data[n == n_plot & delta > 0 & censoring_rate == 0.5],
  truth = 3,
  breaks = seq(-0.05, 10.5, by = 0.25),
  facet = facet_thesis(rows = "model", cols = "delta", scales = "free_y"),
  x_limits = c(0, 10)
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
  results_dir
)

save_plot(p_type1, "type1_error", plot_dir, width = 8, height = 4.5)
save_plot(p_power, "power_curv_cens50", plot_dir, width = 10, height = 6.5)
save_plot(p_box, "splot_points_boxplot_cens50", plot_dir, width = 10, height = 6.5)
save_plot(p_hist, "split_points_hist_cens50n500", plot_dir, width = 10, height = 5)
