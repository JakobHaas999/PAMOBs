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

# Plot results
dodge <- position_dodge(width = 0.4)

p <- ggplot(
  summary_table,
  aes(
    x = scenario,
    y = rejection_rate,
    colour = factor(n),
    group = factor(n)
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "grey45"
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = dodge,
    width = 0.15
  ) +
  geom_point(
    position = dodge,
    size = 2.8
  ) +
  facet_grid(
    rows = vars(covariate_setting, df),
    cols = vars(censoring_rate),
    labeller = labeller(
      covariate_setting = c(
        none = "Without covariates",
        included = "With covariates"
      ),
      censoring_rate = c(
        "0.3" = "30% censoring",
        "0.5" = "50% censoring",
        "0.7" = "70% censoring"
      )
    )
  ) +
  scale_x_discrete(
    limits = c("cubic", "weibull", "lognormal"),
    labels = c("Cubic", "Weibull", "Lognormal")
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    limits = c(0, NA)
  ) +
  scale_colour_manual(
    values = c("500" = "#0072B2", "1000" = "#D55E00"),
    name = "Sample size"
  ) +
  labs(
    x = "Scenario",
    y = "Split rate"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p)

# Plot cubic and weibull hazard with df = 6
summary_table <- subset(summary_table, (df == 6) & (scenario %in% c("cubic", "weibull")))
