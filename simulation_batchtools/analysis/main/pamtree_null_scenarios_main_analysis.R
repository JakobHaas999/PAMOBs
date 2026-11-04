library("data.table")
library("batchtools")
library("ggplot2")

reg <- suppressWarnings(
  loadRegistry(
    "simulation_batchtools/registries/main/pamtree_null_scenarios_main_v1/",
    work.dir = getwd(),
    writeable = FALSE
  )
)

res <- reduceResultsDataTable()
res <- res[, result[[1]], by = "job.id"]
results <- merge(unwrap(getJobPars()), res, by = "job.id")

summary_h0 <- results[,
  {
    ok <- is.na(error) & is.finite(p_root)
    n_ok <- sum(ok)

    n_splits <- sum(split[ok])
    split_rate <- if (n_ok) n_splits / n_ok else NA_real_
    ci <- if (n_ok) binom.test(n_splits, n_ok)$conf.int else c(NA_real_, NA_real_)

    list(
      n_planned = .N,
      n_completed = n_ok,
      split_rate = split_rate,
      ci_lower = ci[1],
      ci_upper = ci[2]
    )
  },
  by = c("variant", "censoring_rate")
]
print(summary_h0)

# Plot Type-I-error rate
summary_h0[, variant := factor(variant, levels = c("oracle_log", "sp0", "sp30", "REML"))]
dodge <- position_dodge(width = 0.35)
p_type1 <- ggplot(
  summary_h0,
  aes(
    x = factor(censoring_rate),
    y = split_rate,
    colour = variant,
    group = variant
  )
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
  geom_point(position = dodge, size = 3) +
  scale_colour_manual(
    values = c(oracle_log = "#0072B2", sp0 = "#009E73", sp30 = "#E69F00", REML = "#D55E00"),
    labels = c(
      oracle_log = "Oracle: log(t)", sp0 = "P-spline, sp = 0",
      sp30 = "P-spline, sp = 30", REML = "P-spline, REML"
    )
  ) +
  scale_colour_manual(
    values = c(oracle_log = "#000000", sp0 = "#0072B2", sp30 = "#D55E00", REML = "#CC79A7"),
    labels = c(
      oracle_log = "Oracle: log(t)", sp0 = "P-spline, sp = 0",
      sp30 = "P-spline, sp = 30", REML = "P-spline, REML"
    )
  ) +
  scale_x_discrete(labels = function(x) scales::percent(as.numeric(x), accuracy = 1)) +
  scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  labs(
    title = "Rejection rate under H0 when splitting on time",
    subtitle = "Weibull baseline hazard | n = 1000",
    x = "Target censoring rate",
    y = "Rejection Rate",
    colour = "Nodemodel"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )
print(p_type1)

### Possible further analyses (not done) -------------------------------
# - ECDF of root p-values per node model (one censoring rate):
#   shows the same as the rejection rates, but over the whole p-value range.
# - Histogram of tau_hat for false splits per node model:
#   only about one third of false splits lie at t < 2, so no clear
#   concentration at early times (where log(t) is steepest).
# - Distribution of REML-selected sp in the root node (sp_root):
#   median ~180, 10%-90% quantile ~48-544, i.e. much stiffer than sp = 30.
#   Reported in the text instead of a figure.

### Save tables and plots ----------------------------------------------
saveRDS(results, file.path(results_dir, "job_results.rds"))
fwrite(summary_h0, file.path(results_dir, "summary_h0.csv"))

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
