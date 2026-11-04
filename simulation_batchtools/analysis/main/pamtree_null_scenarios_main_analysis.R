### Analysis: node-model diagnosis under H0 (split on time) -------------
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "pamtree_null_scenarios_main_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

# One row per job and node model (variant)
results <- load_results(run_id)
mark_valid(results)

### Checks
sum(!is.na(results$error))
censoring_check <- check_censoring(results[variant == "oracle_log"], by = "censoring_rate")
censoring_check

### Type-I-Error rate per node model
summary_h0 <- summarise_rejection(
  results[valid == TRUE],
  by = c("variant", "censoring_rate")
)
print(summary_h0)

p_type1 <- plot_rejection_by_censoring(
  summary_h0,
  group = "variant",
  y_max = 1,
  title = "Rejection rate under H0 when splitting on time",
  subtitle = "Weibull baseline hazard | n = 1000"
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
save_tables(list(censoring_check = censoring_check, summary_h0 = summary_h0), results_dir)
save_plot(p_type1, "type1_error", plot_dir)
