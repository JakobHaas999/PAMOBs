### Analysis: PEMtree null scenarios (split on time, H0) ----------------
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "null_scenario_main_df_v2"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

# Registry was deleted; results are read from the cached job_results.rds
metrics <- load_results(run_id)
metrics[, covariate_setting := factor(covariate_setting, levels = c("none", "included"))]

issues <- metrics[!(success %in% TRUE), .(
  job.id, covariate_setting, scenario, rep, seed, completed, error
)]
if (nrow(issues)) {
  print(issues)
  warning("Not all jobs completed successfully", call. = FALSE)
}

design_cols <- c("n", "censoring_rate", "covariate_setting", "scenario", "df")
summary_table <- summarise_rejection(
  metrics[success %in% TRUE],
  by = design_cols,
  split_col = "split_detected"
)
summary_table <- merge(
  summary_table,
  check_censoring(metrics[success %in% TRUE], by = design_cols)[, .(
    n, censoring_rate, covariate_setting, scenario, df, mean_observed_censoring = observed
  )],
  by = design_cols
)
print(summary_table)
save_tables(list(summary = summary_table), results_dir)

# Only discuss cubic and weibull hazard for now
summary_table <- summary_table[scenario %in% c("cubic", "weibull")]

## Split rates ---------------------------------------------------
rejection_plots <- lapply(c("weibull", "cubic"), function(s) {
  plot_rejection_by_censoring(
    summary_table[scenario == s],
    group = "n",
    facet = facet_thesis(rows = "covariate_setting", cols = "df"),
    title = labels_scenario[[s]]
  )
})
names(rejection_plots) <- c("weibull", "cubic")
print(rejection_plots$weibull)
print(rejection_plots$cubic)

## Root p values -------------------------------------------
n_plot <- 1000
censoring_plot <- 0.5

p_ecdf <- plot_pvalue_ecdf(
  metrics[success %in% TRUE & is.finite(root_p_value) &
    scenario %in% c("cubic", "weibull") & n == n_plot & censoring_rate == censoring_plot],
  group = "df",
  p_col = "root_p_value",
  facet = facet_thesis(rows = "covariate_setting", cols = "scenario"),
  title = "Distribution of root-test p-values",
  subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%")
)
p_ecdf

## Locations of false splits ---------------------------------
tau_base <- metrics[
  scenario %in% c("cubic", "weibull") & n == n_plot & censoring_rate == censoring_plot
]
tau_base[, df := factor(df, levels = c(8, 6, 4))]
tau_base[, false_split := success %in% TRUE & split_detected %in% TRUE & is.finite(tau_hat)]

tau_data <- tau_base[false_split == TRUE]
tau_counts <- tau_base[, .(k = sum(false_split)), by = .(scenario, covariate_setting, df)]

p_tau <- ggplot(tau_data, aes(x = tau_hat, y = df)) +
  geom_point(
    position = position_jitter(width = 0, height = 0.14, seed = 1),
    colour = "grey25", alpha = 0.65, size = 2
  ) +
  geom_text(
    data = tau_counts,
    aes(x = 9.8, label = paste0("k = ", k)),
    position = position_nudge(y = 0.28),
    hjust = 1, colour = "grey40", size = 3
  ) +
  facet_thesis(rows = "covariate_setting", cols = "scenario", drop = FALSE) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  scale_y_discrete(labels = label_df, drop = FALSE) +
  labs(
    title = "Locations of false splits under the null hypothesis",
    subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%"),
    x = expression("Estimated split time " * hat(tau)),
    y = NULL
  ) +
  theme_thesis(base_size = 12)
p_tau

# Save plots ---------------------------------------------------------------
save_plot(rejection_plots$weibull, "H0_rejection_rates_weibull", plot_dir, width = 11, height = 6)
save_plot(rejection_plots$cubic, "H0_rejection_rates_cubic", plot_dir, width = 11, height = 6)
save_plot(p_ecdf, "H0_root_pvalues_ecdf", plot_dir, width = 10, height = 7)
save_plot(p_tau, "H0_false_split_locations", plot_dir, width = 10, height = 7)
