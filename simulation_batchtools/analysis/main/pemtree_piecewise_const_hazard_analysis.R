### Analysis: Piecewise constant hazard (split on time, correctly specified node model)
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "pemtree_piecewise_constant_hazard_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")
tau <- 5

## Load results
res <- load_results(run_id)
mark_valid(res)

## Checks
invalid_jobs <- res[, .(n_invalid = sum(valid == FALSE, share_invalid = mean(valid == FALSE))),
  by = "delta"
][n_invalid > 0]
invalid_jobs
censoring_check <- check_censoring(res, by = c("delta", "censoring_rate"))
censoring_check

## Type I error and power
summary_rejection <- summarise_rejection(res[valid == TRUE], by = "delta")
summary_rejection

p_ecdf <- plot_pvalue_ecdf(
  res[delta == 0 & valid == TRUE],
  group = "n",
  subtitle = "delta = 0 | n = 500 | 50% censoring"
)
p_ecdf

## Split points (true tau = 5)
split_data <- res[valid == TRUE & split == TRUE]
summary_split_point <- summarise_split_points(
  split_data[delta > 0],
  truth = tau, by = "delta", within = 0.5
)
summary_split_point

p_hist <- plot_split_points_hist(
  split_data[delta > 0],
  truth = tau, breaks = seq(-0.05, 10.05, by = 0.1),
  facet = facet_thesis_wrap("delta", ncol = 1, scales = "free_y"),
  fill = pal_n[["500"]], x_breaks = 0:10, x_lab = expression(hat(tau)),
  subtitle = "n = 500 | 50% censoring | dashed line: true break"
)
p_hist

## Save
save_tables(list(
  invalid_jobs = invalid_jobs, censoring_check = censoring_check,
  summary_rejection = summary_rejection, summary_split_point = summary_split_point
), results_dir)

save_plot(p_ecdf, "pvalues_ecdf_h0", plot_dir)
save_plot(p_hist, "split_points_hist", plot_dir, width = 9, height = 6)
