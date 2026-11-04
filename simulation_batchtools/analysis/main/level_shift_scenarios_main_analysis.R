### Analysis: PEMtree level shift in time (split on time) ---------------
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "level_shift_main_df_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

# Read from the cached job_results.rds (the registry is 3.3 GB)
metrics <- load_results(run_id)
metrics[, covariate_setting := factor(covariate_setting, levels = c("none", "included"))]

## Issues
issues <- metrics[!(success %in% TRUE), .(
  job.id, scenario, n, df, delta, censoring_rate,
  covariate_setting, rep, seed, completed, error
)]
if (nrow(issues)) {
  print(issues)
  warning("Not all jobs completed successfully.", call. = FALSE)
}

missing_root_tests <- metrics[
  success %in% TRUE & !is.finite(root_p_value),
  .(
    job.id, scenario, n, df, delta, censoring_rate,
    covariate_setting, rep, seed, completed, root_p_value, split_detected
  )
]

## Create summary table
design_cols <- c("scenario", "n", "df", "delta", "censoring_rate", "covariate_setting")
summary_rejection <- summarise_rejection(
  metrics[success %in% TRUE],
  by = design_cols,
  split_col = "split_detected"
)
# Location accuracy is conditional on a successful fit with a detected split (true tau = 3)
summary_tau <- summarise_split_points(
  metrics[success %in% TRUE & split_detected %in% TRUE & is.finite(tau_hat)],
  truth = 3,
  by = design_cols,
  split_col = "tau_hat"
)
summary_table <- merge(summary_rejection, summary_tau, by = design_cols, all.x = TRUE, suffixes = c("", "_tau"))
print(summary_table)

save_tables(
  list(issues = issues, missing_root_tests = missing_root_tests, summary = summary_table),
  results_dir
)

## Plots

# Power with confidence interval for different bs df, one plot per hazard ratio
hrs <- c(1.25, 1.5, 2, 3)
power_plots <- list()
for (s in c("weibull", "cubic")) {
  for (hr in hrs) {
    power_plots[[paste0("power_", s, "_hr_", gsub("\\.", "_", hr))]] <- plot_rejection_by_censoring(
      summary_table[abs(delta - log(hr)) < 1e-8 & scenario == s],
      group = "n",
      facet = facet_thesis(rows = "covariate_setting", cols = "df"),
      y_max = 1,
      ref_line = NULL,
      title = labels_scenario[[s]],
      subtitle = sprintf("delta = %.3f | Hazard ratio = %.2f", log(hr), hr)
    )
  }
}
power_plots$power_weibull_hr_3

# Boxplots of estimated change points (df = 6, without covariates)
tau_data <- metrics[
  success %in% TRUE & split_detected %in% TRUE & is.finite(tau_hat) &
    df == 6 & covariate_setting == "none"
]
box_plots <- lapply(c(1.5, 3), function(hr) {
  plot_split_points_box(
    tau_data[abs(delta - log(hr)) < 1e-8],
    truth = 3,
    x = "censoring_rate",
    split_col = "tau_hat",
    facet = facet_thesis_wrap("scenario"),
    x_lab = "Target censoring rate",
    x_labels = label_percent_x,
    y_lab = expression(hat(tau)),
    subtitle = paste0("Hazard ratio = ", hr, " | df = 6 | Covariate setting: none"),
    title = NULL
  ) +
    labs(caption = "Successful fits with a detected split only. Dashed line: true change point.") +
    theme(plot.caption = element_text(hjust = 0))
})
names(box_plots) <- paste0("split_locations_boxplot_hr_", gsub("\\.", "_", c(1.5, 3)), "_df_6_none")
box_plots[[2]]

# Histogram of estimated change points (HR = 3, df = 6, none, n = 1000, 30 % censoring)
p_tau_hist <- plot_split_points_hist(
  tau_data[abs(delta - log(3)) < 1e-8 & n == 1000 & abs(censoring_rate - 0.3) < 1e-8],
  truth = 3,
  breaks = seq(-0.05, 10.05, by = 0.1),
  split_col = "tau_hat",
  facet = facet_thesis_wrap("scenario"),
  x_lab = expression(hat(tau)),
  subtitle = "Hazard ratio = 3 | df = 6 | Covariate setting: none\nn = 1000 | Target censoring rate = 30%"
) +
  labs(
    y = "Number of detected splits",
    caption = "Successful fits with a detected split only. Dashed line: true change point."
  ) +
  theme(plot.caption = element_text(hjust = 0))
p_tau_hist

## Save plots
for (nm in names(power_plots)) save_plot(power_plots[[nm]], nm, plot_dir, width = 10, height = 7)
for (nm in names(box_plots)) save_plot(box_plots[[nm]], nm, plot_dir)
save_plot(p_tau_hist, "split_locations_histogram_hr_3_df_6_none_n_1000_censor_0_3", plot_dir)
