### setup
source("setup.R")
registry_dir <- "simulation_batchtools/registries/main/time_changing_effect_main_v1"
results_dir <- "simulation_batchtools/results/main/time_changing_effect_main_v1"

reg <- batchtools::loadRegistry(
  file.dir = registry_dir, work.dir = getwd(), conf.file = NA,
  writeable = FALSE, make.default = FALSE
)
job_design <- readRDS(file.path(reg$file.dir, "job_design.rds"))

done_ids <- batchtools::findDone(reg = reg)
collected <- batchtools::reduceResultsDataTable(
  ids = done_ids, reg = reg
)

metrics <- rbindlist(collected$result, fill = TRUE)
metrics[, job.id := collected$job.id]
metrics <- merge(as.data.table(job_design), metrics, by = "job.id", all.x = TRUE)
metrics[, completed := job.id %in% collected$job.id]

issues <- metrics[
  !completed | success %in% FALSE,
  .(job.id, scenario, seed, completed, error)
]

if (nrow(issues)) {
  print(issues)
  warning("Not all jobs completed successfully.", call. = FALSE)
}

# Save individual results and issues before creating plots
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(metrics, file.path(results_dir, "job_results.rds"))
data.table::fwrite(issues, file.path(results_dir, "issues.csv"))

# Create a summary table
scenario_cols <- c("scenario", "delta", "censoring_rate", "n", "df")
## 1) Null scenario
summary_h0 <- metrics[delta == 0,
  {
    ok <- success %in% TRUE
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    split_rate <- if (n_ok) n_splits / n_ok else NA_real_
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    list(
      n_planned = .N,
      n_success = n_ok,
      n_splits = n_splits,
      split_rate = split_rate,
      ci_lower = ci[1],
      ci_upper = ci[2]
    )
  },
  by = scenario_cols
]
print(summary_h0)
data.table::fwrite(summary_h0, file.path(results_dir, "summary_h0.csv"))

# Plots
dodge <- position_dodge(width = 0.35)
p_h0 <- ggplot(
  summary_h0[n != 250],
  aes(
    x = factor(censoring_rate),
    y = split_rate,
    colour = factor(n),
    shape = factor(n),
    group = factor(n)
  )
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
  geom_point(position = dodge, size = 3) +
  facet_wrap(~scenario, labeller = labeller(
    scenario = c(
      "cubic" = "Cubic baseline hazard",
      "weibull" = "Weibull baseline hazard"
    )
  )) +
  scale_colour_manual(values = c("500" = "#009E73", "1000" = "#7B3294")) +
  scale_shape_manual(values = c("500" = 16, "1000" = 17)) +
  scale_x_discrete(
    labels = function(x) scales::percent(as.numeric(x), accuracy = 1)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    x = "Target censoring rate",
    y = "Rejection Rate",
    colour = "Sample size",
    shape = "Sample size"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

print(p_h0)

# Plot root p value
n_plot <- 1000
censoring_plot <- 0.4
p_data <- metrics[
  success %in% TRUE & is.finite(root_p_value) & delta == 0 &
    n == n_plot & censoring_rate == censoring_plot
]
p_data$scenario <- factor(p_data$scenario, levels = c("cubic", "weibull"))

p_ecdf <- ggplot(p_data, aes(x = root_p_value)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dashed") +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  facet_wrap(
    ~scenario,
    labeller = labeller(
      scenario = c(
        "cubic" = "Cubic baseline hazard",
        "weibull" = "Weibull baseline hazard"
      )
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Distribution of root-test p-values",
    subtitle = paste0(
      "n = ", n_plot,
      " | Target censoring rate = ", 100 * censoring_plot, "%"
    ),
    x = "Root p-value",
    y = "Cumulative proportion"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, colour = "grey40")
  )

print(p_ecdf)

# 2) Time-changing covariate effect
summary_table_shift <- metrics[delta != 0,
  {
    ok <- success %in% TRUE
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    power <- if (n_ok) n_splits / n_ok else NA_real_
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    mean_tau_hat <- if (n_splits > 0) mean(tau_hat, na.rm = TRUE) else NA_real_
    mean_abs_tau_error <- if (n_splits) mean(abs_tau_error, na.rm = TRUE) else NA_real_

    list(
      n_planned = .N,
      n_success = n_ok,
      n_splits = n_splits,
      power = power,
      ci_lower = ci[1],
      ci_upper = ci[2],
      mean_tau_hat = mean_tau_hat,
      mean_abs_tau_error = mean_abs_tau_error
    )
  },
  by = scenario_cols
]
print(summary_table_shift)
data.table::fwrite(summary_table_shift, file.path(results_dir, "summary_power.csv"))

# Plot empirical power
p_power <- ggplot(
  summary_table_shift[n != 250],
  aes(
    x = factor(censoring_rate),
    y = power,
    colour = factor(n),
    shape = factor(n),
    group = factor(n)
  )
) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = dodge, width = 0.1, linewidth = 0.6
  ) +
  geom_point(position = dodge, size = 3) +
  facet_grid(
    rows = vars(delta),
    cols = vars(scenario),
    labeller = labeller(
      delta = as_labeller(
        c(
          "0.3" = "delta == 0.3",
          "0.6" = "delta == 0.6",
          "0.9" = "delta == 0.9"
        ),
        default = label_parsed
      ),
      scenario = c(
        cubic = "Cubic baseline hazard",
        weibull = "Weibull baseline hazard"
      )
    )
  ) +
  labs(
    x = "Target censoring rate",
    y = "Power",
    colour = "Sample size",
    shape = "Sample size"
  ) +
  theme_bw() +
  scale_colour_manual(
    values = c("500" = "#009E73", "1000" = "#7B3294")
  ) +
  scale_shape_manual(
    values = c("500" = 16, "1000" = 17)
  ) +
  scale_x_discrete(
    labels = function(x) scales::percent(as.numeric(x), accuracy = 1)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = grid::unit(1.2, "lines")
  )

print(p_power)

# Analysis of split points
# Split point Analysis

# Boxplots
tau_data <- metrics[
  success %in% TRUE & split_detected %in% TRUE & is.finite(tau_hat) &
    delta > 0 & n %in% c(500, 1000) & df == 6
]

p_tau <- ggplot(
  tau_data,
  aes(
    x = factor(censoring_rate),
    y = tau_hat,
    fill = factor(n),
    group = interaction(censoring_rate, n)
  )
) +
  geom_hline(
    yintercept = 3,
    linetype = "dashed",
    colour = "grey45"
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    width = 0.55,
    colour = "grey30",
    outlier.alpha = 0.25,
    outlier.size = 1
  ) +
  facet_grid(
    rows = vars(delta),
    cols = vars(scenario),
    labeller = labeller(
      delta = as_labeller(
        c(
          "0.3" = "delta == 0.3",
          "0.6" = "delta == 0.6",
          "0.9" = "delta == 0.9"
        ),
        default = label_parsed
      ),
      scenario = c(
        cubic = "Cubic baseline hazard",
        weibull = "Weibull baseline hazard"
      )
    )
  ) +
  scale_fill_manual(values = c("500" = "#009E73", "1000" = "#7B3294")) +
  scale_x_discrete(
    labels = function(x) scales::percent(as.numeric(x), accuracy = 1)
  ) +
  labs(
    x = "Target censoring rate",
    y = expression(hat(tau)),
    fill = "Sample size",
    caption = paste(
      "Successful fits with a detected split only.",
      "Dashed line: true change point at t = 3."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

print(p_tau)

# Histogram
delta_plot <- 0.6
tau_data_hist <- tau_data[
  censoring_rate == censoring_plot & n == n_plot & delta == delta_plot
]
p_tau_hist <- ggplot(
  tau_data_hist,
  aes(
    x = tau_hat
  )
) +
  geom_histogram(
    binwidth = 0.1,
    center = 3,
    fill = "#0072B2",
    colour = "white",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 3, linetype = "dashed", colour = "red"
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
  labs(
    x = expression(hat(tau)),
    y = "Number of detected splits",
    subtitle = paste0(
      "n = ", n_plot,
      " | Target censoring rate = ", 100 * censoring_plot, "%",
      " | delta = ", delta_plot
    )
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_tau_hist)

# Save plots
plot_dir <- file.path(results_dir, "figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, filename, width = 8, height = 5) {
  # Vector graphic for LaTeX
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )

  # High-resolution PNG
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

save_plot(p_h0, "type1_error", width = 8, height = 5)
save_plot(p_ecdf, "pvalue_ecdf", width = 8, height = 5)
save_plot(p_power, "power", width = 9, height = 8)
save_plot(p_tau, "split_locations_boxplot", width = 9, height = 8)
save_plot(p_tau_hist, "split_locations_histogram", width = 8, height = 5)
