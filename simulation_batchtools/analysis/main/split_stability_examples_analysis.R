### : Analysis split stability examples
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "split_stability_examples_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

res <- load_results(run_id)

# Summary table
stability <- res[,
  {
    sp <- split_point[split %in% TRUE]
    .(
      split_rate = mean(split, na.rm = TRUE), k = length(sp),
      full_point = full_split_point[1], median = median(sp),
      lower = quantile(sp, 0.025), upper = quantile(sp, 0.975)
    )
  },
  by = c("setting", "delta", "censoring_rate")
]
stability[, width := upper - lower]
stability[, reliability := fifelse(k < 250, "Fewer than 250 splits", "At least 250 splits")]

save_tables(list(summary_stability = stability), results_dir)

# Plot Intervals
plot_stability_intervals <- function(data, truth, y_lab, y_limits, y_breaks = waiver()) {
  ggplot(data, aes(x = factor(censoring_rate), colour = reliability)) +
    geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.8) +
    geom_point(aes(y = median, shape = "Subsample median"), size = 2.5) +
    geom_point(
      aes(y = full_point, shape = "Split on the full data"),
      position = position_nudge(x = 0.18), size = 2.5, stroke = 0.9, na.rm = TRUE
    ) +
    geom_hline(yintercept = truth, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
    geom_text(
      data = data[k < 250], aes(y = upper, label = paste0("k = ", k)),
      nudge_y = diff(y_limits) * 0.04, size = 3, colour = "grey30", show.legend = FALSE
    ) +
    facet_thesis_wrap("delta") +
    scale_colour_manual(
      values = c("At least 250 splits" = pal_n[["1000"]], "Fewer than 250 splits" = "grey60"),
      name = NULL
    ) +
    scale_shape_manual(values = c("Subsample median" = 16, "Split on the full data" = 4), name = NULL) +
    scale_x_discrete(labels = label_percent_x) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    labs(
      x = "Censoring", y = y_lab,
      subtitle = "Lines: 2.5% to 97.5% percentiles over subsamples with a split | dashed line: true split point"
    ) +
    theme_thesis(discrete_x = TRUE)
}

z_limits <- c(-0.5, 0.75)
z_breaks <- seq(-0.5, 0.75, by = 0.25)

save_plot(
  plot_stability_intervals(stability[setting == "covariate"], 0, "Split point in z", z_limits, z_breaks),
  "stability_interval_covariate",
  plot_dir,
  width = 8,
  height = 4.5
)
save_plot(
  plot_stability_intervals(stability[setting == "time_changing_effect"], 3, "Split point in t", c(0, 6), 0:6),
  "stability_interval_time_changing_effect",
  plot_dir,
  width = 8,
  height = 4.5
)
save_plot(
  plot_stability_intervals(stability[setting == "structural_change"], 3, "Split point in t", c(0, 6), 0:6),
  "stability_interval_structural_change",
  plot_dir,
  width = 8,
  height = 4
)

plot_stability_hist <- function(cell_setting, cell_delta, cell_censoring, truth, breaks,
                                x_lab, x_limits, x_breaks) {
  cell <- res[setting == cell_setting & delta == cell_delta & censoring_rate == cell_censoring]
  plot_split_points_hist(
    cell[split == TRUE],
    truth = truth, breaks = breaks,
    x_lab = x_lab, x_limits = x_limits, x_breaks = x_breaks,
    subtitle = sprintf(
      "delta = %s, %d%% censoring | %d of %d subsamples with a split | dashed line: true split point",
      cell_delta, 100 * cell_censoring, sum(cell$split, na.rm = TRUE), nrow(cell)
    )
  )
}

t_breaks <- seq(-0.05, 10.05, by = 0.1)
z_bins <- seq(-1.025, 1.025, by = 0.05)
save_plot(
  plot_stability_hist("time_changing_effect", 1, 0.5, 3, t_breaks, "Split point in t", c(0, 6), 0:6),
  "stability_hist_time_changing_effect", plot_dir,
  width = 7, height = 3.5
)
save_plot(
  plot_stability_hist("structural_change", 1, 0.5, 3, t_breaks, "Split point in t", c(0, 6), 0:6),
  "stability_hist_structural_change", plot_dir,
  width = 7, height = 3.5
)
save_plot(
  plot_stability_hist("covariate", 0.5, 0.5, 0, z_bins, "Split point in z", z_limits, z_breaks),
  "stability_hist_covariate", plot_dir,
  width = 7, height = 3.5
)
