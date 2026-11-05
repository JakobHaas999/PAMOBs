### Analysis: Structural change in the hazard (split on time)
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "structural_change_hazard_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

res <- load_results(run_id)
mark_valid(res)

## 1. Checks
invalid_jobs <- res[, .(n_valid = sum(valid == FALSE), share_invalid = mean(valid == FALSE)),
  by = c("scenario", "df", "n", "censoring_rate")
][n_valid > 0]
censoring_check <- check_censoring(res[df == 6], by = c("scenario", "delta", "censoring_rate"))

## 2. Type I error (delta = 0)
summary_h0 <- summarise_rejection(res[delta == 0 & valid == TRUE],
  by = c("scenario", "df", "n", "censoring_rate")
)
p_type1 <- plot_rejection_by_censoring(
  summary_h0,
  group = "n",
  facet = facet_thesis(rows = "scenario", cols = "df")
)
p_type1

p_ecdf <- plot_pvalue_ecdf(
  res[delta == 0 & valid == TRUE & n == 1000 & censoring_rate == 0.5],
  group = "df", facet = facet_thesis_wrap("scenario"),
  subtitle = "n = 1000 | 50% censoring"
)
p_ecdf

## 3. Power
summary_power <- summarise_rejection(
  res[valid == TRUE],
  by = c("scenario", "df", "n", "censoring_rate", "delta")
)
summary_power

power_curves <- lapply(sort(unique(summary_power$censoring_rate)), function(rate) {
  plot_power_curve(
    summary_power[censoring_rate == rate],
    facet = facet_thesis(rows = "scenario", cols = "df"),
    x_breaks = unique(summary_power$delta),
    subtitle = paste0(100 * rate, "% censoring | HR = exp(delta): 1, 1.65, 2.72, 4.48")
  )
})
names(power_curves) <- paste0("power_curves_cens", c(30, 50, 70))

## 4. Split points (true tau = 3)
split_data <- res[valid == TRUE & split == TRUE & delta > 0]
summary_split_point <- summarise_split_points(
  split_data,
  truth = 3,
  by = c("scenario", "df", "n", "censoring_rate", "delta"),
  split_col = "tau_hat"
)

p_box <- plot_split_points_box(
  split_data[censoring_rate == 0.3],
  truth = 3, split_col = "tau_hat",
  facet = facet_thesis(rows = "scenario", cols = "df"),
  y_limits = c(0, 10), y_breaks = seq(0, 10, by = 2),
  y_lab = expression(hat(tau)),
  subtitle = "30% censoring | only replications with a split | dashed line: true break"
)
p_box

p_hist <- plot_split_points_hist(
  split_data[delta == 1.5 & n == 1000 & censoring_rate == 0.3],
  truth = 3, breaks = seq(-0.05, 10.05, by = 0.1), split_col = "tau_hat",
  facet = facet_thesis(rows = "df", cols = "scenario", scales = "free_y"),
  x_lab = expression(hat(tau)),
  subtitle = "delta = 1.5 | n = 1000 | 30% censoring"
)

## 5. Save
save_tables(list(
  invalid_jobs = invalid_jobs, censoring_check = censoring_check,
  summary_h0 = summary_h0, summary_power = summary_power,
  summary_split_point = summary_split_point
), results_dir)

save_plot(p_type1, "type1_error", plot_dir, width = 10, height = 6.5)
save_plot(p_ecdf, "pvalues_ecdf_n1000_cens50", plot_dir)
for (nm in names(power_curves)) save_plot(power_curves[[nm]], nm, plot_dir, width = 10, height = 6.5)
save_plot(p_box, "split_points_box_cens30", plot_dir, width = 10, height = 6.5)
save_plot(p_hist, "split_points_hist_delta1.5_n1000_cens30", plot_dir, width = 10, height = 8)
