##### Analysis of 02 level shift null

source("setup.R")

run_id <- "level_shift_null_robust_v1"
results_folder <- file.path("simulation/results/main", run_id, "blocks")
blocks <- list.files(results_folder, pattern = "^block_[0-9]+\\.rds$")
results <- rbindlist(lapply(blocks, function(b) {
  as.data.table(readRDS(file.path(results_folder, b)))
}))

size_summary <- results[, {
  n_success <- sum(success)
  n_splits <- sum(split_detected[success])
  test <- binom.test(n_splits, n_success, p = 0.05)

  list(
    n_runs = .N,
    n_errors = sum(!n_success),
    n_splits = n_splits,
    alpha_hat = n_splits / n_success,
    ci_low = test$conf.int[1],
    ci_high = test$conf.int[2],
    p_value = test$p.value
  )
}, by = c("n", "censoring_rate")]

size_summary