### Anaylsis of PEM tree with a lognormal baseline and comparison with a weibull AFT tree
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "pemtree_covariate_lognormal_v2"
results_dir <- results_path(run_id)
plots_dir <- file.path(results_dir, "figures")

res <- load_results(run_id)
mark_valid(res)

## 1. Checks
invalid_jobs <- res[, .(n_invalid = sum(!valid), share_invalid = mean(!valid)),
  by = c("n", "delta", "censoring_rate", "model")
][n_invalid > 0]
invalid_jobs
censoring_check <- check_censoring(res, by = c("n", "delta", "censoring_rate", "model"))
censoring_check

## 2. Type I error
summary_h0 <- summarise_rejection(
  res[delta == 0 & valid == TRUE],
  by = c("n", "delta", "censoring_rate", "model")
)
summary_h0

p_type1 <- plot_rejection_by_censoring(
  summary_h0,
  group = "n",
  facet = facet_thesis(rows = "model")
)
p_type1
