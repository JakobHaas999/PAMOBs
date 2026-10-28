library("ggplot2")

registry_dir <- "simulation_batchtools/registries/main/tree_comparison_main_v1"
results_dir <- "simulation_batchtools/results/main/tree_comparison_main_v1"

reg <- batchtools::loadRegistry(
  file.dir = registry_dir,
  conf.file = NA,
  writeable = FALSE,
  make.default = FALSE
)

job_design <- readRDS(file.path(reg$file.dir, "job_design.rds"))

collected <- batchtools::reduceResultsDataTable(
  ids = batchtools::findDone(reg = reg),
  reg = reg
)

metrics <- data.table::rbindlist(collected$result, fill = TRUE)
metrics[, job.id := collected$job.id]
metrics <- merge(data.table::as.data.table(job_design), metrics, by = "job.id")

summary_table <- metrics[,
  {
    ok <- success %in% TRUE &
      test_available %in% TRUE &
      is.finite(root_p_value) &
      !is.na(split_detected)

    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    list(
      n_completed = .N,
      n_valid = n_ok,
      n_failed = sum(success %in% FALSE),
      n_missing_test = sum(
        success %in% TRUE & !(test_available %in% TRUE)
      ),
      n_splits = n_splits,
      split_rate = if (n_ok) n_splits / n_ok else NA_real_,
      ci_lower = ci[[1]],
      ci_upper = ci[[2]]
    )
  },
  by = .(scenario, censoring_rate, delta, n, df, method)
]

# Save results
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(metrics, file.path(results_dir, "job_results.rds"))
data.table::fwrite(summary_table, file.path(results_dir, "summary.csv"))

comparison <- data.table::dcast(
  summary_table,
  scenario + censoring_rate + delta + n + df ~ method,
  value.var = "split_rate"
)
print(comparison)

# Plots
# 1) Null scenario

summary_table_0 <- summary_table[delta == 0]
dodge <- position_dodge(width = 0.35)
p0 <- ggplot(
  summary_table_0,
  aes(
    x = factor(censoring_rate),
    y = split_rate,
    colour = factor(method),
    shape = factor(method),
    group = factor(method)
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = dodge,
    width = 0.10,
    linewidth = 0.6
  ) +
  geom_point(
    position = dodge,
    size = 3
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
  scale_colour_manual(
    values = c(pemtree = "#0072B2", PErMtree = "#D55E00")
  ) +
  scale_shape_manual(
    values = c(pemtree = 16, PErMtree = 17)
  ) +
  scale_x_discrete(
    labels = function(x) scales::percent(as.numeric(x), accuracy = 1)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = "Target censoring rate",
    y = "Rejection rate",
    colour = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, colour = "grey40"),
    panel.spacing = grid::unit(1.2, "lines")
  )

print(p0)

# 2) H1 scenario
summary_table_shift <- summary_table[delta != 0]
p_shift <- ggplot(
  summary_table_shift,
  aes(
    x = factor(censoring_rate), y = split_rate,
    colour = factor(method), shape = factor(method), group = factor(method)
  )
) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
  geom_point(position = dodge, size = 3) +
  facet_grid(
    rows = vars(delta),
    cols = vars(scenario),
    labeller = labeller(
      delta = c(
        "0.4" = "delta = 0.4",
        "0.8" = "delta = 0.8"
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
    colour = "Method",
    shape = "Method"
  ) +
  theme_bw() +
  scale_colour_manual(
    values = c(pemtree = "#0072B2", PErMtree = "#D55E00")
  ) +
  scale_x_discrete(
    labels = function(x) scales::percent(as.numeric(x), accuracy = 1)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    breaks = seq(0.75, 1, by = 0.05),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, colour = "grey40"),
    panel.spacing = grid::unit(1.2, "lines")
  )

print(p_shift)

# Plot distribution of p values under H0
p_data <- subset(metrics, success %in% TRUE & is.finite(root_p_value) & delta == 0)
p_data$censoring_rate <- factor(p_data$censoring_rate)

p_ecdf <- ggplot(
  p_data,
  aes(x = root_p_value, colour = method, group = method)
) +
  geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dashed") +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  facet_wrap(
    ~censoring_rate,
    labeller = labeller(
      censoring_rate = c(
        "0.4" = "Target censoring rate = 0.4",
        "0.6" = "target censoring rate = 0.6"
      )
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.25)
  ) +
  scale_y_continuous(
    limits = c(0, 1), labels = scales::label_percent(accuracy = 1)
  ) +
  scale_colour_manual(
    values = c(pemtree = "#0072B2", PErMtree = "#D55E00")
  ) +
  labs(
    title = "Distribution of root-test p-values",
    x = "Root p-value",
    y = "Cumulative proportion",
    group = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p_ecdf)

# Save plots
plot_dir <- file.path(results_dir, "figures")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot, filename, width = 8, height = 5) {
  # Vector graphic for LaTeX
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )

  # High-resolution PNG
  ggsave(
    filename = file.path(plot_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

save_plot(p0, "type1_error", width = 8, height = 5)
save_plot(p_shift, "power", width = 9, height = 6)
save_plot(p_ecdf, "pvalue_ecdf", width = 8, height = 5)
