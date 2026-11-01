source("setup.R")
registry_dir <- "simulation_batchtools/registries/main/level_shift_main_df_v1"
results_dir <- "simulation_batchtools/results/main/level_shift_main_df_v1"

reg <- batchtools::loadRegistry(
  file.dir = registry_dir,
  work.dir = getwd(),
  conf.file = NA,
  make.default = FALSE,
  writeable = FALSE
)
job_design <- readRDS(file.path(reg$file.dir, "job_design.rds"))

metrics <- data.table::data.table(
  job.id = integer(), success = logical(), error = character(),
  split_detected = logical(), n_terminal_nodes = integer(),
  tau_hat = numeric(), abs_tau_error = numeric(),
  observed_censoring = numeric(), root_p_value = numeric()
)
done_ids <- batchtools::findDone(reg = reg)

if (nrow(done_ids)) {
  collected <- batchtools::reduceResultsDataTable(ids = done_ids, reg = reg)
  collected_metrics <- data.table::rbindlist(collected$result, fill = TRUE)
  collected_metrics[, job.id := collected$job.id]
  metrics <- data.table::rbindlist(list(metrics, collected_metrics), fill = TRUE)
}

metrics <- merge(
  job_design, metrics,
  by = "job.id", all.x = TRUE
)
metrics <- data.table::as.data.table(metrics)
metrics[, completed := job.id %in% done_ids$job.id]

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
  .(job.id, scenario, n, df, delta, censoring_rate,
    covariate_setting, rep, seed, completed, root_p_value, split_detected)
]

# Save all individual results before creating summaries and plots
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(metrics, file.path(results_dir, "job_results.rds"))
data.table::fwrite(issues, file.path(results_dir, "issues.csv"))
data.table::fwrite(missing_root_tests, file.path(results_dir, "missing_root_tests.csv"))

## Create summary table
scenario_cols <- setdiff(colnames(job_design), c("job.id", "rep", "seed"))
summary_table <- metrics[,
  {
    ok <- success %in% TRUE
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    # Location accuracy is conditional on a successful fit with a detected split
    located <- ok & split_detected %in% TRUE & is.finite(tau_hat)
    tau_error <- tau_hat[located] - tau[1L]
    mean_tau_hat <- if (any(located)) mean(tau_hat[located]) else NA_real_
    mean_abs_tau_error <- if (any(located)) mean(abs(tau_error)) else NA_real_

    list(
      n_planned = .N,
      n_completed = sum(completed),
      n_success = n_ok,
      n_failed = sum(success %in% FALSE),
      n_missing = sum(!completed),
      n_missing_test = sum(ok & !is.finite(root_p_value)),
      n_splits = n_splits,
      split_rate = if (n_ok > 0) n_splits / n_ok else NA_real_,
      ci_lower = ci[1],
      ci_upper = ci[2],
      n_tau = sum(located),
      mean_tau_hat = mean_tau_hat,
      mean_abs_tau_error = mean_abs_tau_error,
      bias_tau = if (any(located)) mean(tau_error) else NA_real_,
      rmse_tau = if (any(located)) sqrt(mean(tau_error^2)) else NA_real_,
      mean_observed_censoring = if (n_ok) mean(observed_censoring[ok]) else NA_real_
    )
  },
  by = scenario_cols
]
print(summary_table)
data.table::fwrite(summary_table, file.path(results_dir, "summary.csv"))

## Plots

# Plot power with confidence interval for different bs df
plot_power <- function(summary_table, delta_plot, scenario_plot) {
  plot_data <- data.table::as.data.table(summary_table)[
    abs(delta - delta_plot) < 1e-8 &
      scenario == scenario_plot
  ]

  if (!nrow(plot_data)) {
    stop("No results available for the selected settings.")
  }

  plot_data[, covariate_setting := factor(
    covariate_setting,
    levels = c("none", "included")
  )]

  dodge <- position_dodge(width = 0.35)

  ggplot(
    plot_data,
    aes(
      x = factor(censoring_rate),
      y = split_rate,
      colour = factor(n),
      shape = factor(n),
      group = factor(n)
    )
  ) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = dodge,
      width = 0.1,
      linewidth = 0.6
    ) +
    geom_point(position = dodge, size = 3) +
    facet_grid(
      rows = vars(covariate_setting),
      cols = vars(df),
      labeller = labeller(
        df = function(x) paste0("df = ", x),
        covariate_setting = c(
          none = "Without covariates",
          included = "With covariates"
        )
      )
    ) +
    scale_colour_manual(
      values = c(
        "250" = "#E69F00",
        "500" = "#009E73",
        "1000" = "#7B3294"
      )
    ) +
    scale_shape_manual(
      values = c("250" = 15, "500" = 16, "1000" = 17)
    ) +
    scale_x_discrete(
      labels = function(x) {
        scales::percent(as.numeric(x), accuracy = 1)
      }
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    labs(
      title = switch(scenario_plot,
        cubic = "Cubic baseline hazard",
        weibull = "Weibull baseline hazard",
        scenario_plot
      ),
      subtitle = sprintf(
        "delta = %.3f | Hazard ratio = %.2f",
        delta_plot, exp(delta_plot)
      ),
      x = "Target censoring rate",
      y = "Power",
      colour = "Sample size",
      shape = "Sample size"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = grid::unit(1.2, "lines")
    )
}

p0_weibull <- plot_power(summary_table, delta_plot = log(1.25), scenario_plot = "weibull")
print(p0_weibull)
p1_weibull <- plot_power(summary_table, delta_plot = log(1.5), scenario_plot = "weibull")
print(p1_weibull)
p2_weibull <- plot_power(summary_table, delta_plot = log(2), scenario_plot = "weibull")
print(p2_weibull)
p3_weibull <- plot_power(summary_table, delta_plot = log(3), scenario_plot = "weibull")
print(p3_weibull)

p0_cubic <- plot_power(summary_table, delta_plot = log(1.25), scenario_plot = "cubic")
print(p0_cubic)
p1_cubic <- plot_power(summary_table, delta_plot = log(1.5), scenario_plot = "cubic")
print(p1_cubic)
p2_cubic <- plot_power(summary_table, delta_plot = log(2), scenario_plot = "cubic")
print(p2_cubic)
p3_cubic <- plot_power(summary_table, delta_plot = log(3), scenario_plot = "cubic")
print(p3_cubic)

# Plot boxplot of taus
plot_tau_boxplot <- function(
  metrics,
  delta_plot,
  df_plot = 6,
  covariate_setting_plot = "none"
) {
  plot_data <- data.table::as.data.table(metrics)[
    success %in% TRUE &
      split_detected %in% TRUE &
      is.finite(tau_hat) &
      abs(delta - delta_plot) < 1e-8 &
      df == df_plot &
      covariate_setting == covariate_setting_plot
  ]

  if (!nrow(plot_data)) {
    stop("No detected splits available for the selected settings.")
  }

  true_tau <- unique(plot_data[, .(scenario, tau)])

  ggplot(
    plot_data,
    aes(
      x = factor(censoring_rate),
      y = tau_hat,
      fill = factor(n),
      group = interaction(censoring_rate, n)
    )
  ) +
    geom_hline(
      data = true_tau,
      aes(yintercept = tau),
      inherit.aes = FALSE,
      linetype = "dashed",
      colour = "grey45"
    ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.65,
      colour = "grey30",
      outlier.alpha = 0.25,
      outlier.size = 1
    ) +
    facet_wrap(
      ~scenario,
      labeller = labeller(
        scenario = c(
          cubic = "Cubic baseline hazard",
          weibull = "Weibull baseline hazard"
        )
      )
    ) +
    scale_fill_manual(
      values = c(
        "250" = "#E69F00",
        "500" = "#009E73",
        "1000" = "#7B3294"
      )
    ) +
    scale_x_discrete(
      labels = function(x) {
        scales::percent(as.numeric(x), accuracy = 1)
      }
    ) +
    labs(
      subtitle = paste0(
        "Hazard ratio = ", round(exp(delta_plot), 2),
        " | df = ", df_plot,
        " | Covariate setting: ", covariate_setting_plot
      ),
      x = "Target censoring rate",
      y = expression(hat(tau)),
      fill = "Sample size",
      caption = paste(
        "Successful fits with a detected split only.",
        "Dashed line: true change point."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0)
    )
}

p1_box <- plot_tau_boxplot(metrics, delta_plot = log(1.5))
print(p1_box)
p2_box <- plot_tau_boxplot(metrics, delta_plot = log(3))
print(p2_box)


plot_tau_histogram <- function(
  metrics,
  delta_plot,
  df_plot = 6,
  covariate_setting_plot = "none",
  n_plot = 1000,
  censoring_plot = 0.3,
  binwidth = 0.1
) {
  plot_data <- data.table::as.data.table(metrics)[
    success %in% TRUE &
      split_detected %in% TRUE &
      is.finite(tau_hat) &
      abs(delta - delta_plot) < 1e-8 &
      df == df_plot &
      covariate_setting == covariate_setting_plot &
      n == n_plot &
      abs(censoring_rate - censoring_plot) < 1e-8
  ]

  if (!nrow(plot_data)) {
    stop("No detected splits available for the selected settings.")
  }

  true_tau <- unique(plot_data[, .(scenario, tau)])

  ggplot(
    plot_data,
    aes(x = tau_hat, fill = factor(n))
  ) +
    geom_histogram(
      binwidth = binwidth,
      center = 0,
      colour = "white",
      linewidth = 0.2,
      show.legend = FALSE
    ) +
    geom_vline(
      data = true_tau,
      aes(xintercept = tau),
      linetype = "dashed",
      colour = "red",
      linewidth = 0.6
    ) +
    facet_wrap(
      ~scenario,
      labeller = labeller(
        scenario = c(
          cubic = "Cubic baseline hazard",
          weibull = "Weibull baseline hazard"
        )
      )
    ) +
    scale_fill_manual(
      values = c(
        "250" = "#E69F00",
        "500" = "#009E73",
        "1000" = "#7B3294"
      )
    ) +
    labs(
      subtitle = paste0(
        "Hazard ratio = ", round(exp(delta_plot), 2),
        " | df = ", df_plot,
        " | Covariate setting: ", covariate_setting_plot,
        "\nn = ", n_plot,
        " | Target censoring rate = ", 100 * censoring_plot, "%"
      ),
      x = expression(hat(tau)),
      y = "Number of detected splits",
      caption = paste(
        "Successful fits with a detected split only.",
        "Dashed line: true change point."
      )
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0)
    )
}

p_tau_hist <- plot_tau_histogram(metrics, delta_plot = log(3))
print(p_tau_hist)

## Save plots
level_shift_plot_dir <- file.path(results_dir, "figures")
dir.create(level_shift_plot_dir, recursive = TRUE, showWarnings = FALSE)

save_level_shift_plot <- function(plot, filename, width = 9, height = 6) {
  # Vector graphic for LaTeX
  ggsave(
    filename = file.path(level_shift_plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )

  # High-resolution PNG
  ggsave(
    filename = file.path(level_shift_plot_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

save_level_shift_plot(p0_weibull, "power_weibull_hr_1_25", width = 10, height = 7)
save_level_shift_plot(p1_weibull, "power_weibull_hr_1_5", width = 10, height = 7)
save_level_shift_plot(p2_weibull, "power_weibull_hr_2", width = 10, height = 7)
save_level_shift_plot(p3_weibull, "power_weibull_hr_3", width = 10, height = 7)
save_level_shift_plot(p0_cubic, "power_cubic_hr_1_25", width = 10, height = 7)
save_level_shift_plot(p1_cubic, "power_cubic_hr_1_5", width = 10, height = 7)
save_level_shift_plot(p2_cubic, "power_cubic_hr_2", width = 10, height = 7)
save_level_shift_plot(p3_cubic, "power_cubic_hr_3", width = 10, height = 7)
save_level_shift_plot(p1_box, "split_locations_boxplot_hr_1_5_df_6_none")
save_level_shift_plot(p2_box, "split_locations_boxplot_hr_3_df_6_none")
save_level_shift_plot(
  p_tau_hist,
  "split_locations_histogram_hr_3_df_6_none_n_1000_censor_0_3"
)
