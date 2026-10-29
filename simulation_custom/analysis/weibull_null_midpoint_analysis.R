source("setup.R")
### ANalysis of results
results_dir <- "simulation_custom/results"
results <- readRDS(file.path(results_dir, "weibull_null_midpoint.rds"))

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
print(summary)

# Plot the results
model_order <- c(
  "reference_time",
  "reference_mid_time",
  "reference_time_cut",
  "reference_time_all"
)

model_labels <- c(
  reference_time = "Interval end\nStep = 0.1",
  reference_mid_time = "Midpoint\nStep = 0.1",
  reference_time_cut = "Interval end\nStep = 0.5",
  reference_time_all = "Interval end\nEvent-time grid"
)

summary[, model := factor(model, levels = model_order)]

n_rep <- data.table::uniqueN(results$rep)

p_null <- ggplot(
  summary,
  aes(x = model, y = rejection_rate)
) +
  geom_hline(
    yintercept = 0.05,
    colour = "grey55",
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.12,
    linewidth = 0.6,
    colour = "#0072B2"
  ) +
  geom_point(
    size = 3,
    colour = "#0072B2"
  ) +
  scale_x_discrete(labels = model_labels) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1)
  ) +
  coord_cartesian(
    ylim = c(0, max(0.10, summary$ci_upper, na.rm = TRUE))
  ) +
  labs(
    title = "Time discretization and splits under the Weibull null scenario",
    subtitle = paste0(
      "n = 1,000 | ", n_rep,
      " recorded repetitions | Target censoring = 40%"
    ),
    x = "Time evaluation and interval grid",
    y = "Split rate"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )

p_null

# Save the plot
figure_dir <- "simulation_custom/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(
    figure_dir,
    "weibull_null_midpoint_comparison.pdf"
  ),
  plot = p_null,
  width = 9,
  height = 5,
  units = "in"
)
