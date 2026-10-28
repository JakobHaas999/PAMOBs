source("setup.R")

registry_dir <- "simulation_batchtools/registries/main/null_scenario_main_df_v2"
results_dir <- "simulation_batchtools/results/main/null_scenario_main_df_v2"

reg <- batchtools::loadRegistry(
  file.dir = registry_dir, writeable = FALSE, conf.file = NA, make.default = FALSE
)
job_design <- readRDS(file.path(reg$file.dir, "job_design.rds"))

metrics <- data.table::data.table(
  job.id = integer(), success = logical(), error = character(),
  split_detected = logical(), n_terminal_nodes = integer(),
  tau_hat = numeric(), abs_tau_error = numeric(),
  observed_censoring = numeric()
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

issues <- metrics[!(success %in% TRUE), .(
  job.id, covariate_setting, scenario, rep, seed, completed, error
)]
if (nrow(issues)) {
  print(issues)
  warning("Non all jobs completed successfully", call. = FALSE)
}

scenario_cols <- setdiff(colnames(job_design), c("job.id", "rep", "seed"))
summary_table <- metrics[,
  {
    ok <- success %in% TRUE
    n_ok <- sum(ok)
    n_splits <- sum(split_detected[ok])

    rejection_rate <- if (n_ok) n_splits / n_ok else NA_real_

    ci <- if (n_ok) binom.test(n_splits, n_ok, p = 0.05)$conf.int else c(NA_real_, NA_real_)

    list(
      n_planned = .N,
      n_completed = sum(completed),
      n_success = sum(ok),
      n_failed = sum(completed & !ok),
      n_missing = sum(!completed),
      n_splits = n_splits,
      rejection_rate = rejection_rate,
      ci_lower = ci[1],
      ci_upper = ci[2],
      mean_observed_censoring = if (n_ok) mean(observed_censoring[ok]) else NA_real_
    )
  },
  by = scenario_cols
]

summary_table <- summary_table[, .(
  n, censoring_rate, covariate_setting, scenario, df, n_planned,
  n_success, n_failed, n_missing, rejection_rate, ci_lower,
  ci_upper, mean_observed_censoring
)]
print(summary_table)

# Save results
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(metrics, file.path(results_dir, "job_results.rds"))
data.table::fwrite(summary_table, file.path(results_dir, "summary.csv"))

# Only discuss cubic and weibull hazard for now
summary_table <- summary_table[scenario %in% c("cubic", "weibull")]
summary_table

## Split rates ---------------------------------------------------
# 1) Weibull
summary_table_weibull <- summary_table[scenario == "weibull"]
dodge <- position_dodge(width = 0.35)
p_weibull <- ggplot(
  summary_table_weibull,
  aes(
    x = factor(censoring_rate),
    y = rejection_rate,
    colour = factor(n),
    shape = factor(n),
    group = factor(n)
  )
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
  geom_point(position = dodge, size = 3) +
  facet_grid(
    rows = vars(covariate_setting),
    cols = vars(df),
    labeller = labeller(
      covariate_setting = c(
        none = "Without covariates",
        included = "With covariates"
      ),
      df = c(
        "4" = "df = 4",
        "6" = "df = 6",
        "8" = "df = 8"
      )
    )
  ) +
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
    shape = "Sample size",
    title = "Weibull baseline hazard"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

print(p_weibull)

# 2) Cubic
summary_table_cubic <- summary_table[scenario == "cubic"]
dodge <- position_dodge(width = 0.35)
p_cubic <- ggplot(
  summary_table_cubic,
  aes(
    x = factor(censoring_rate),
    y = rejection_rate,
    colour = factor(n),
    shape = factor(n),
    group = factor(n)
  )
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
  geom_point(position = dodge, size = 3) +
  facet_grid(
    rows = vars(covariate_setting),
    cols = vars(df),
    labeller = labeller(
      covariate_setting = c(
        none = "Without covariates",
        included = "With covariates"
      ),
      df = c(
        "4" = "df = 4",
        "6" = "df = 6",
        "8" = "df = 8"
      )
    )
  ) +
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
    shape = "Sample size",
    title = "Cubic baseline hazard"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold")
  )

print(p_cubic)

## Root p values -------------------------------------------
n_plot <- 1000
censoring_plot <- 0.5

p_data <- subset(
  metrics,
  success %in% TRUE &
    is.finite(root_p_value) &
    scenario %in% c("cubic", "weibull") &
    n == n_plot &
    censoring_rate == censoring_plot
)

p_data$covariate_setting <- factor(
  p_data$covariate_setting,
  levels = c("none", "included")
)

p_ecdf <- ggplot(
  p_data,
  aes(x = root_p_value, linetype = factor(df))
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    colour = "grey60",
    linetype = "dashed"
  ) +
  stat_ecdf(
    geom = "step",
    linewidth = 0.8
  ) +
  facet_grid(
    rows = vars(covariate_setting),
    cols = vars(scenario),
    labeller = labeller(
      scenario = c(
        cubic = "Cubic baseline hazard",
        weibull = "Weibull baseline hazard"
      ),
      covariate_setting = c(
        none = "Without covariates",
        included = "With covariates"
      )
    )
  ) +
  scale_linetype_manual(
    values = c("4" = "solid", "6" = "longdash", "8" = "dotdash")
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
    y = "Cumulative proportion",
    linetype = "Spline df"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, colour = "grey40")
  )

p_ecdf

# Plot wrong split points
n_plot <- 1000
censoring_plot <- 0.5

tau_base <- metrics[
  scenario %in% c("cubic", "weibull") &
    n == n_plot &
    censoring_rate == censoring_plot
]

tau_base[, covariate_setting := factor(
  covariate_setting,
  levels = c("none", "included")
)]
tau_base[, df := factor(df, levels = c(8, 6, 4))]

tau_base[, false_split :=
  success %in% TRUE &
    split_detected %in% TRUE &
    is.finite(tau_hat)]

tau_data <- tau_base[false_split == TRUE]

tau_counts <- tau_base[
  , .(k = sum(false_split)),
  by = .(scenario, covariate_setting, df)
]

p_tau <- ggplot(tau_data, aes(x = tau_hat, y = df)) +
  geom_point(
    position = position_jitter(width = 0, height = 0.14, seed = 1),
    colour = "grey25",
    alpha = 0.65,
    size = 2
  ) +
  geom_text(
    data = tau_counts,
    aes(x = 9.8, label = paste0("k = ", k)),
    position = position_nudge(y = 0.28),
    hjust = 1,
    colour = "grey40",
    size = 3
  ) +
  facet_grid(
    rows = vars(covariate_setting),
    cols = vars(scenario),
    drop = FALSE,
    labeller = labeller(
      scenario = c(
        cubic = "Cubic baseline hazard",
        weibull = "Weibull baseline hazard"
      ),
      covariate_setting = c(
        none = "Without covariates",
        included = "With covariates"
      )
    )
  ) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, by = 2)
  ) +
  scale_y_discrete(
    labels = function(x) paste0("df = ", x),
    drop = FALSE
  ) +
  labs(
    title = "Locations of false splits under the null hypothesis",
    subtitle = paste0(
      "n = ", n_plot,
      " | Target censoring rate = ", 100 * censoring_plot, "%"
    ),
    x = expression("Estimated split time " * hat(tau)),
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, colour = "grey40")
  )

p_tau

# Save plots ---------------------------------------------------------------

plot_dir <- file.path(results_dir, "figures")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plot, filename, width = 9, height = 6) {
  # PDF for LaTeX
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

# 1) Rejection rates under H0
save_plot(
  p_weibull,
  "H0_rejection_rates_weibull",
  width = 11,
  height = 6
)

save_plot(
  p_cubic,
  "H0_rejection_rates_cubic",
  width = 11,
  height = 6
)

# 2) Distribution of root-test p-values
save_plot(
  p_ecdf,
  "H0_root_pvalues_ecdf",
  width = 10,
  height = 7
)

# 3) Locations of false splits
save_plot(
  p_tau,
  "H0_false_split_locations",
  width = 10,
  height = 7
)
