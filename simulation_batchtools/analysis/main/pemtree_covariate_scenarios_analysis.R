### Analysis of results ----------------------------------------------
### setup
library("batchtools")
library("data.table")
library("ggplot2")

results_dir <- "simulation_batchtools/results/main/pemtree_covariate_main_v1"
cache <- file.path(results_dir, "job_results.rds")
if (file.exists(cache)) {
  results <- readRDS(cache)
} else {
  reg <- suppressWarnings(
    loadRegistry(
      "simulation_batchtools/registries/main/pemtree_covariate_main_v1",
      work.dir = getwd(),
      writeable = FALSE
    )
  )
  results <- merge(
    unwrap(getJobPars()),
    unwrap(reduceResultsDataTable()),
    by = "job.id"
  )
  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, cache)
}
results[, p_root := vapply(p_root, function(x) if (length(x) == 1) x else NA_real_, numeric(1))]
results[, valid := is.na(error) & !is.na(p_root)]
results[valid == FALSE, split := NA]

### Checks
# Report that
sum(!is.na(results$error))
nrow(results[valid == FALSE])
invalid_jobs <- results[valid == FALSE, .(job.id, scenario, split_type, n, delta, censoring_rate, n_events)]
invalid_jobs
# Observed censoring
censoring_check <- results[, .(
  observed = mean(1 - n_events / n)
), by = c("scenario", "split_type", "delta", "censoring_rate")][
  ,
  diff := observed - censoring_rate
][]
censoring_check

### Type-I-Error rate
summary_0 <- results[delta == 0 & valid == TRUE,
  {
    n_splits <- sum(split)
    split_rate <- n_splits / .N
    ci <- binom.test(n_splits, .N)$conf.int
    list(
      n_valid = .N,
      split_rate = split_rate,
      ci_lower = ci[1],
      ci_upper = ci[2]
    )
  },
  by = c("scenario", "split_type", "n", "censoring_rate")
]
print(summary_0)

# Plot Type-I-error
dodge <- position_dodge(width = 0.35)
p_type1 <- ggplot(
  summary_0,
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
  facet_grid(
    rows = vars(split_type),
    cols = vars(scenario),
    labeller = labeller(
      split_type = c(binary = "Binary z", numeric = "Numeric z"),
      scenario = c(cubic = "Cubic baseline", weibull = "Weibull baseline")
    )
  ) +
  scale_colour_manual(values = c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")) +
  scale_shape_manual(values = c("250" = 15, "500" = 16, "1000" = 17)) +
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

print(p_type1)

# Plot p value ecdf
censoring_plot <- 0.5

p_data <- results[
  delta == 0 & is.na(error) & is.finite(p_root) & censoring_rate == censoring_plot
]
p_ecdf <- ggplot(
  p_data,
  aes(x = p_root, colour = factor(n), linetype = factor(n))
) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dashed") +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  facet_grid(
    rows = vars(split_type),
    cols = vars(scenario),
    labeller = labeller(
      split_type = c(binary = "Binary z", numeric = "Numeric z"),
      scenario = c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard")
    )
  ) +
  scale_colour_manual(values = c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")) +
  scale_linetype_manual(values = c("250" = "solid", "500" = "longdash", "1000" = "dotdash")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "Distribution of root-test p-values under H0",
    subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "%"),
    x = "Root p-value",
    y = "Cumulative proportion",
    colour = "Sample size",
    linetype = "Sample size"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )
print(p_ecdf)

### Power analysis
summary_power <- results[delta > 0 & valid == TRUE,
  {
    n_splits <- sum(split)
    ci <- binom.test(n_splits, .N)$conf.int
    list(n_valid = .N, split_rate = n_splits / .N, ci_lower = ci[1], ci_upper = ci[2])
  },
  by = c("scenario", "split_type", "n", "censoring_rate", "delta")
]
print(summary_power)

# Function to plot power
plot_power <- function(summary_table, delta_plot) {
  plot_data <- summary_table[abs(delta - delta_plot) < 1e-8]
  dodge <- position_dodge(width = 0.35)

  ggplot(
    plot_data,
    aes(
      x = factor(censoring_rate), y = split_rate,
      colour = factor(n), shape = factor(n), group = factor(n)
    )
  ) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
    geom_point(position = dodge, size = 3) +
    facet_grid(
      rows = vars(split_type),
      cols = vars(scenario),
      labeller = labeller(
        split_type = c(binary = "Binary z", numeric = "Numeric z"),
        scenario = c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard")
      )
    ) +
    scale_colour_manual(values = c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")) +
    scale_shape_manual(values = c("250" = 15, "500" = 16, "1000" = 17)) +
    scale_x_discrete(labels = function(x) scales::percent(as.numeric(x), accuracy = 1)) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    labs(
      title = "Power to detect effect modification of x by z",
      subtitle = sprintf(
        "delta = %.2f | Hazard ratio for x: %.2f (z <= 0) vs. %.2f (z > 0)",
        delta_plot, exp(0.6), exp(0.6 + delta_plot)
      ),
      x = "Target censoring rate", y = "Rejection Rate",
      colour = "Sample size", shape = "Sample size"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold")
    )
}

power_plots <- lapply(c(0.25, 0.5, 0.75, 1), function(d) plot_power(summary_power, d))
names(power_plots) <- paste0("power_delta_", c(0.25, 0.5, 0.75, 1))
power_plots$power_delta_0.5

### Split point
tau_hat_data <- results[split_type == "numeric" & valid & split == TRUE]
summary_split_point <- tau_hat_data[, .(
  n_splits = .N,
  mean_split = mean(split_point),
  sd_split = sd(split_point),
  mean_abs_error = mean(abs(split_point)),
  share_within_0.1 = mean(abs(split_point) <= 0.1)
), by = c("scenario", "n", "censoring_rate", "delta")]
summary_split_point

# Plot split points as boxplot
censoring_plot <- 0.5

p_box <- ggplot(
  tau_hat_data[censoring_rate == censoring_plot],
  aes(x = factor(delta), y = split_point, fill = factor(n))
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(
    cols = vars(scenario),
    labeller = labeller(scenario = c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard"))
  ) +
  scale_fill_manual(values = c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
  labs(
    title = "Estimated split points (numeric z, true split at 0)",
    subtitle = paste0("Target censoring rate = ", 100 * censoring_plot, "% | only replications with a split"),
    x = expression(delta), y = "Estimated split point", fill = "Sample size"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

p_box

# Histogram
n_plot <- 1000

p_hist <- ggplot(
  tau_hat_data[n == n_plot & censoring_rate == censoring_plot],
  aes(x = split_point)
) +
  geom_histogram(breaks = seq(-1, 1, by = 0.05), fill = "#7B3294", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  facet_grid(
    rows = vars(delta),
    cols = vars(scenario),
    labeller = labeller(
      delta = function(x) paste0("delta = ", x),
      scenario = c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard")
    ),
    scales = "free_y"
  ) +
  labs(
    title = "Distribution of estimated split points (numeric z)",
    subtitle = paste0("n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot, "%"),
    x = "Estimated split point", y = "Count"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )
p_hist

### Save tables and plots ----------------------------------------------
fwrite(invalid_jobs, file.path(results_dir, "invalid_jobs.csv"))
fwrite(censoring_check, file.path(results_dir, "censoring_check.csv"))
fwrite(summary_0, file.path(results_dir, "summary_h0.csv"))
fwrite(summary_power, file.path(results_dir, "summary_power.csv"))
fwrite(summary_split_point, file.path(results_dir, "summary_split_point.csv"))

plot_dir <- file.path(results_dir, "figures")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot, filename, width = 9, height = 6) {
  # PDF for LaTeX
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "in", device = cairo_pdf
  )
  # High-resolution PNG
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "in", dpi = 300
  )
}

save_plot(p_type1, "type1_error")
save_plot(p_ecdf, paste0("p_value_ecdf_cens", 100 * censoring_plot))
for (nm in names(power_plots)) save_plot(power_plots[[nm]], nm)
save_plot(p_box, paste0("split_points_boxplot_cens", 100 * censoring_plot))
save_plot(p_hist, paste0("split_points_hist_n", n_plot, "_cens", 100 * censoring_plot), height = 9)
