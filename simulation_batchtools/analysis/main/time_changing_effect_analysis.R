library("batchtools")
library("data.table")
library("ggplot2")

results_dir <- "simulation_batchtools/results/main/time_changing_effect_main_v2"
cache <- file.path(results_dir, "job_results.rds")
if (file.exists(cache)) {
  results <- readRDS(cache)
} else {
  reg <- suppressWarnings(
    loadRegistry(
      "simulation_batchtools/registries/main/time_changing_effect_main_v2",
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


### Checks
sum(!is.na(results$error))
results[, valid := is.na(error) & !is.na(p_root)]
results[valid == FALSE, split := NA]

invalid_jobs <- results[valid == FALSE, .(job.id, scenario, n, delta, censoring_rate, n_events)]
invalid_jobs

# Observed censoring
censoring_check <- results[, .(
  observed = mean(1 - n_events / n)
), by = c("scenario", "delta", "n", "censoring_rate")]
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
  by = c("scenario", "n", "censoring_rate")
][
  order(scenario, n, censoring_rate)
]
summary_0

# Plot Type-I-Error
dodge <- position_dodge(width = 0.35)
p_type1 <- ggplot(
  summary_0,
  aes(
    x = factor(censoring_rate),
    y = split_rate,
    colour = factor(n),
    group = factor(n)
  )
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
  geom_point(position = dodge, size = 3) +
  facet_wrap(
    ~scenario,
    labeller = labeller(
      scenario = c(cubic = "Cubic baseline", weibull = "Weibull baseline")
    )
  ) +
  scale_color_manual(values = c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")) +
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
    colour = "Sample size"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

p_type1


### Power analysis
summary_power <- results[delta > 0 & valid == TRUE,
  {
    n_splits <- sum(split)
    ci <- binom.test(n_splits, .N)$conf.int
    list(n_valid = .N, split_rate = n_splits / .N, ci_lower = ci[1], ci_upper = ci[2])
  },
  by = c("scenario", "n", "censoring_rate", "delta")
][
  order(scenario, n, censoring_rate, delta)
]
summary_power

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
    facet_wrap(
      ~scenario,
      labeller = labeller(
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
      title = "Power to detect a time-varying effect of x",
      subtitle = sprintf(
        "delta = %.2f | Hazard ratio for x: %.2f (t < 3) vs. %.2f (t >= 3)",
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
power_plots$power_delta_0.75
power_plots$power_delta_1

### Power curves over delta
power_curve_data <- rbind(
  copy(summary_0)[, delta := 0],
  summary_power,
  use.names = TRUE
)

p_power_curve <- ggplot(
  power_curve_data,
  aes(x = delta, y = split_rate, colour = factor(n), shape = factor(n), group = factor(n))
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_line(linewidth = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.03, linewidth = 0.5) +
  geom_point(size = 2.5) +
  facet_grid(
    rows = vars(scenario),
    cols = vars(censoring_rate),
    labeller = labeller(
      scenario = c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard"),
      censoring_rate = function(x) paste0(100 * as.numeric(x), "% censoring")
    )
  ) +
  scale_colour_manual(values = c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")) +
  scale_shape_manual(values = c("250" = 15, "500" = 16, "1000" = 17)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1)
  ) +
  labs(
    title = "Power to detect a time-varying effect of x",
    subtitle = "Effect of x: 0.6 before t = 3 and 0.6 + delta afterwards",
    x = expression(delta), y = "Rejection Rate",
    colour = "Sample size", shape = "Sample size"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

p_power_curve


### Split point analysis
tau_hat_data <- results[valid == TRUE & split == TRUE & delta > 0]
summary_split_point <- tau_hat_data[, .(
  n_splits = .N,
  mean_split = mean(split_point),
  sd_split = sd(split_point),
  mean_abs_error = mean(abs(split_point - tau))
), by = c("scenario", "n", "censoring_rate", "delta")][
  order(scenario, n, censoring_rate, delta)
]
summary_split_point

# Plot split points as bocplot
censoring_plot <- 0.5

p_box <- ggplot(
  tau_hat_data[censoring_rate == censoring_plot],
  aes(x = factor(delta), y = split_point, fill = factor(n))
) +
  geom_hline(yintercept = 3, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(
    ~scenario,
    labeller = labeller(scenario = c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard"))
  ) +
  scale_fill_manual(values = c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2)) +
  labs(
    title = "Estimated split points in time (true change at t = 3)",
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

# Histogram of split points (Weibull, n = 1000)
n_plot <- 1000

p_hist <- ggplot(
  tau_hat_data[scenario == "weibull" & n == n_plot & censoring_rate == censoring_plot],
  aes(x = split_point)
) +
  geom_histogram(
    breaks = seq(-0.05, 10.05, by = 0.1),
    fill = "#7B3294", colour = "white", linewidth = 0.2
  ) +
  geom_vline(xintercept = 3, linetype = "dashed", colour = "grey45") +
  facet_wrap(
    ~delta,
    ncol = 1,
    scales = "free_y",
    labeller = labeller(delta = function(x) paste0("delta = ", x))
  ) +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
  labs(
    title = "Distribution of estimated split points (Weibull baseline hazard)",
    subtitle = paste0(
      "n = ", n_plot, " | Target censoring rate = ", 100 * censoring_plot,
      "% | true change at t = 3"
    ),
    x = "Estimated split point", y = "Count"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

p_hist



### Coefficients of x before and after the split
coef_data <- results[valid & split == TRUE & delta > 0]
summary_coef <- coef_data[, .(
  n_splits = .N,
  mean_before = mean(x_before),
  bias_before = mean(x_before - beta),
  mean_after = mean(x_after),
  true_after = beta[1] + delta[1],
  bias_after = mean(x_after - (beta + delta)),
  rmse_before = sqrt(mean((x_before - beta)^2)),
  rmse_after = sqrt(mean((x_after - (beta + delta))^2))
), by = c("scenario", "n", "censoring_rate", "delta")][
  order(scenario, n, censoring_rate, delta)
]
summary_coef

# Boxplot of estimates vs. true values (n = 1000, 50 % censoring)
coef_long <- melt(
  coef_data[n == 1000 & censoring_rate == 0.5],
  id.vars = c("scenario", "delta", "beta"),
  measure.vars = c("x_before", "x_after"),
  variable.name = "period", value.name = "estimate"
)
coef_long[, truth := fifelse(period == "x_before", beta, beta + delta)]
coef_long[, period := factor(period,
  levels = c("x_before", "x_after"),
  labels = c("Before t = 3", "After t = 3")
)]
truth_data <- unique(coef_long[, .(scenario, delta, period, truth)])

p_coef <- ggplot(coef_long, aes(x = factor(delta), y = estimate, fill = period)) +
  geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(
    data = truth_data, aes(x = factor(delta), y = truth, group = period),
    position = position_dodge(width = 0.8), shape = 4, size = 3, stroke = 1.2,
    colour = "black", inherit.aes = FALSE
  ) +
  facet_wrap(
    ~scenario,
    labeller = labeller(scenario = c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard"))
  ) +
  scale_fill_manual(values = c("Before t = 3" = "grey80", "After t = 3" = "#56B4E9")) +
  labs(
    title = "Estimated effect of x in the two terminal nodes",
    subtitle = "n = 1000 | Target censoring rate = 50% | only replications with a split | x = true value",
    x = expression(delta), y = "Estimated coefficient of x", fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

p_coef

### Save tables and plots ----------------------------------------------
fwrite(invalid_jobs, file.path(results_dir, "invalid_jobs.csv"))
fwrite(censoring_check, file.path(results_dir, "censoring_check.csv"))
fwrite(summary_0, file.path(results_dir, "summary_h0.csv"))
fwrite(summary_power, file.path(results_dir, "summary_power.csv"))
fwrite(summary_split_point, file.path(results_dir, "summary_split_point.csv"))
fwrite(summary_coef, file.path(results_dir, "summary_coef.csv"))

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
for (nm in names(power_plots)) save_plot(power_plots[[nm]], nm)
save_plot(p_power_curve, "power_curves", width = 10, height = 6.5)
save_plot(p_box, paste0("split_points_boxplot_cens", 100 * censoring_plot))
save_plot(p_hist, paste0("split_points_hist_weibull_n", n_plot, "_cens", 100 * censoring_plot), height = 8)
save_plot(p_coef, "coefficients_x_n1000_cens50")
