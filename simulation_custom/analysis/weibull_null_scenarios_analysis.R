source("setup.R")
## Analysis of results
results_dir <- "simulation_custom/results"
results <- readRDS(file.path(results_dir, "weibull_null_scenarios.rds"))
# Create a summary table
results <- data.table::as.data.table(results)
summary <- results[,
  {
    ok <- success %in% TRUE &
      is.finite(root_p_value)
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    rejection_rate <- if (n_ok) n_splits / n_ok else NA_real_
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    list(
      n_planned = .N, n_success = n_ok, n_split = n_splits,
      rejection_rate = rejection_rate, ci_lower = ci[1], ci_upper = ci[2]
    )
  },
  by = "model"
]

# Plot the results
summary[, baseline_model := factor(
  ifelse(
    grepl("^spline_", model),
    "Spline (df = 6)",
    "Weibull reference"
  ),
  levels = c("Spline (df = 6)", "Weibull reference")
)]
summary[, split_variable := factor(
  ifelse(
    grepl("_time$", model),
    "Time",
    "Independent covariate"
  ),
  levels = c("Time", "Independent covariate")
)]

p_null <- ggplot(
  summary,
  aes(
    x = split_variable,
    y = rejection_rate,
    colour = baseline_model,
    shape = baseline_model,
    group = baseline_model
  )
) +
  geom_hline(
    yintercept = 0.05,
    colour = "grey55",
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = dodge,
    width = 0.12,
    linewidth = 0.6
  ) +
  geom_point(
    position = dodge,
    size = 3
  ) +
  scale_colour_manual(
    values = c(
      "Spline (df = 6)" = "#0072B2",
      "Weibull reference" = "#D55E00"
    )
  ) +
  scale_shape_manual(values = c(16, 17)) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1)
  ) +
  coord_cartesian(
    ylim = c(0, max(0.10, summary$ci_upper, na.rm = TRUE))
  ) +
  labs(
    title = "Splits under the Weibull null scenario",
    subtitle = paste0(
      "n = ", n_obs,
      " | 1000 repetitions",
      " | Target censoring = 40%"
    ),
    x = "Splitting variable",
    y = "Split rate",
    colour = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )
p_null
figure_dir <- "simulation_custom/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(figure_dir, "weibull_null_comparison.pdf"),
  plot = p_null,
  width = 8,
  height = 5,
  units = "in"
)
