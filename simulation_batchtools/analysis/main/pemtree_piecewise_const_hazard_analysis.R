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
  res[delta == 0 & valid == TRUE], group = "n",
  subtitle = "delta = 0 | n = 500 | 50% censoring"
)
p_ecdf

## Split points (true tau = 5)
split_data <- res[valid == TRUE & split ==TRUE]
summary_split_point <- summarise_split_points(
  split_data[delta > 0], truth = tau, by = "delta", within = 0.5
)
summary_split_point


Bash Show template analysis script and helper details
IN
cat -n main/structural_change_hazard_main_analysis.R; sed -n '140,206p;356,385p' analysis_helpers.R; sed -n '78,90p' analysis_helpers.R

OUT
1	### Analysis: Structural change in the hazard (split on time)
2	source("simulation_batchtools/analysis/analysis_helpers.R")
3	
4	run_id <- "structural_change_hazard_v1"
5	results_dir <- results_path(run_id)
6	plot_dir <- file.path(results_dir, "figures")
7	
8	res <- load_results(run_id)
9	mark_valid(res)
10	
11	## 1. Checks
12	invalid_jobs <- res[, .(n_valid = sum(valid == FALSE), share_invalid = mean(valid == FALSE)),
                       13	  by = c("scenario", "df", "n", "censoring_rate")
                       14	][n_valid > 0]
15	censoring_check <- check_censoring(res[df == 6], by = c("scenario", "delta", "censoring_rate"))
16	
17	## 2. Type I error (delta = 0)
18	summary_h0 <- summarise_rejection(res[delta == 0 & valid == TRUE],
                                     19	  by = c("scenario", "df", "n", "censoring_rate")
                                     20	)
21	p_type1 <- plot_rejection_by_censoring(
  22	  summary_h0,
  23	  group = "n",
  24	  facet = facet_thesis(rows = "scenario", cols = "df")
  25	)
26	p_type1
27	
28	p_ecdf <- plot_pvalue_ecdf(
  29	  res[delta == 0 & valid == TRUE & n == 1000 & censoring_rate == 0.5],
  30	  group = "df", facet = facet_thesis_wrap("scenario"),
  31	  subtitle = "n = 1000 | 50% censoring"
  32	)
33	p_ecdf
34	
35	## 3. Power
36	summary_power <- summarise_rejection(
  37	  res[valid == TRUE],
  38	  by = c("scenario", "df", "n", "censoring_rate", "delta")
  39	)
40	summary_power
41	
42	power_curves <- lapply(sort(unique(summary_power$censoring_rate)), function(rate) {
  43	  plot_power_curve(
    44	    summary_power[censoring_rate == rate],
    45	    facet = facet_thesis(rows = "scenario", cols = "df"),
    46	    x_breaks = unique(summary_power$delta),
    47	    subtitle = paste0(100 * rate, "% censoring | HR = exp(delta): 1, 1.65, 2.72, 4.48")
    48	  )
  49	})
50	names(power_curves) <- paste0("power_curves_cens", c(30, 50, 70))
51	
52	## 4. Split points (true tau = 3)
53	split_data <- res[valid == TRUE & split == TRUE & delta > 0]
54	summary_split_point <- summarise_split_points(
  55	  split_data,
  56	  truth = 3,
  57	  by = c("scenario", "df", "n", "censoring_rate", "delta"),
  58	  split_col = "tau_hat"
  59	)
60	
61	p_box <- plot_split_points_box(
  62	  split_data[censoring_rate == 0.3],
  63	  truth = 3, split_col = "tau_hat",
  64	  facet = facet_thesis(rows = "scenario", cols = "df"),
  65	  y_limits = c(0, 10), y_breaks = seq(0, 10, by = 2),
  66	  y_lab = expression(hat(tau)),
  67	  subtitle = "30% censoring | only replications with a split | dashed line: true break"
  68	)
69	p_box
70	
71	p_hist <- plot_split_points_hist(
  72	  split_data[delta == 1.5 & n == 1000 & censoring_rate == 0.3],
  73	  truth = 3, breaks = seq(-0.05, 10.05, by = 0.1), split_col = "tau_hat",
  74	  facet = facet_thesis(rows = "df", cols = "scenario", scales = "free_y"),
  75	  x_lab = expression(hat(tau)),
  76	  subtitle = "delta = 1.5 | n = 1000 | 30% censoring"
  77	)
78	
79	## 5. Save
80	save_tables(list(
  81	  invalid_jobs = invalid_jobs, censoring_check = censoring_check,
  82	  summary_h0 = summary_h0, summary_power = summary_power,
  83	  summary_split_point = summary_split_point
  84	), results_dir)
85	
86	save_plot(p_type1, "type1_error", plot_dir, width = 10, height = 6.5)
87	save_plot(p_ecdf, "pvalues_ecdf_n1000_cens50", plot_dir)
88	for (nm in names(power_curves)) save_plot(power_curves[[nm]], nm, plot_dir, width = 10, height = 6.5)
89	save_plot(p_box, "split_points_box_cens30", plot_dir, width = 10, height = 6.5)
90	save_plot(p_hist, "split_points_hist_delta1.5_n1000_cens30", plot_dir, width = 10, height = 8)
# Mark jobs without a valid root test (error or missing p-value) as invalid.
# Modifies `results` by reference.
mark_valid <- function(results) {
  if (is.list(results$p_root)) {
    results[, p_root := vapply(p_root, function(x) if (length(x) == 1) x else NA_real_, numeric(1))]
  }
  results[, valid := is.na(error) & !is.na(p_root)]
  results[valid == FALSE, split := NA]
  invisible(results)
}


# Summaries ---------------------------------------------------------------
# Observed vs. target censoring rate
check_censoring <- function(results, by) {
  out <- if ("observed_censoring" %in% names(results)) {
    results[, .(observed = mean(observed_censoring)), by = by]
  } else {
    results[, .(observed = mean(1 - n_events / n)), by = by]
  }
  out[, diff := observed - censoring_rate][]
}

# Split (rejection) rate with exact binomial CI; pass only valid jobs.
# Used for the type I error (delta = 0) and for the power (delta > 0).
summarise_rejection <- function(data, by, split_col = "split") {
  out <- data[,
              {
                k <- sum(get(split_col))
                ci <- binom.test(k, .N)$conf.int
                list(n_valid = .N, n_splits = k, split_rate = k / .N, ci_lower = ci[1], ci_upper = ci[2])
              },
              by = by
  ]
  setorderv(out, by)
  out[]
}

# Accuracy of the estimated split point; pass only jobs with a split.
summarise_split_points <- function(data, truth, by, within = NULL, split_col = "split_point") {
  out <- data[,
              {
                sp <- get(split_col)
                res <- list(
                  n_splits = .N,
                  mean_split = mean(sp),
                  sd_split = sd(sp),
                  bias = mean(sp - truth),
                  mean_abs_error = mean(abs(sp - truth)),
                  rmse = sqrt(mean((sp - truth)^2))
                )
                if (!is.null(within)) res[[paste0("share_within_", within)]] <- mean(abs(sp - truth) <= within)
                res
              },
              by = by
  ]
  setorderv(out, by)
  out[]
}


# Plots -------------------------------------------------------------------
# Split rate with CI against the censoring rate.
#   group = "n"       -> colour/shape by sample size
#   group = "variant" -> colour/shape by node model
#   y_max = NULL      -> type I error scale (0 to data max, pretty breaks)
#   y_max = 1         -> power scale (0 to 100 %)
save_plot <- function(plot, filename, dir, width = 9, height = 6) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path(dir, paste0(filename, ".pdf")),
    plot = plot, width = width, height = height, units = "in",
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    filename = file.path(dir, paste0(filename, ".png")),
    plot = plot, width = width, height = height, units = "in",
    dpi = 300, bg = "white"
  )
}

# Save a named list of tables as <name>.csv
save_tables <- function(tables, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(tables)) fwrite(tables[[nm]], file.path(dir, paste0(nm, ".csv")))
}
facet_thesis <- function(rows = NULL, cols = NULL, ...) {
  facet_grid(
    stats::as.formula(paste(if (is.null(rows)) "." else rows, "~", if (is.null(cols)) "." else cols)),
    labeller = thesis_labeller, ...
  )
}

facet_thesis_wrap <- function(var, ...) {
  facet_wrap(stats::as.formula(paste("~", var)), labeller = thesis_labeller, ...)
}

theme_thesis <- function(base_size = 11, discrete_x = FALSE) {
  th <- theme_bw(base_size = base_size) +
    Bash Write and test an analysis skeleton in the scratchpad
  IN
  S=/private/tmp/claude-501/-Users-jakobhaas-Documents-AAA-Bachelorarbeit-PAMOBs/0049965c-7998-43e6-a818-37af68019f2f/scratchpad
  cat > $S/pch_analysis_skeleton.R <<'EOF'
  ### Analysis: Piecewise constant hazard (split on time, correctly specified node model)
  source("simulation_batchtools/analysis/analysis_helpers.R")
  
  run_id <- "pemtree_piecewise_constant_hazard_v1"
  results_dir <- results_path(run_id)
  plot_dir <- file.path(results_dir, "figures")
  tau <- 5
  
  res <- load_results(run_id)
  mark_valid(res)
  
  ## 1. Checks
  invalid_jobs <- res[, .(n_invalid = sum(!valid), share_invalid = mean(!valid)), by = "delta"][n_invalid > 0]
  censoring_check <- check_censoring(res, by = c("delta", "censoring_rate"))
  
  ## 2. Type I error and power (one table, delta = 0 is the type I error)
  summary_rejection <- summarise_rejection(res[valid == TRUE], by = c("n", "censoring_rate", "delta"))
  
  p_ecdf <- plot_pvalue_ecdf(
    res[delta == 0 & valid == TRUE], group = "n",
    subtitle = "delta = 0 | n = 500 | 50% censoring"
  )
  
  ## 3. Split points (true tau = 5)
  split_data <- res[valid == TRUE & split == TRUE]
  summary_split_point <- summarise_split_points(
    split_data[delta > 0], truth = tau, by = c("n", "censoring_rate", "delta"), within = 0.5
  )
  
  p_hist <- plot_split_points_hist(
    split_data[delta > 0], truth = tau, breaks = seq(-0.05, 10.05, by = 0.1),
    facet = facet_thesis_wrap("delta", ncol = 1, scales = "free_y"),
    fill = pal_n[["500"]], x_breaks = 0:10, x_lab = expression(hat(tau)),
    subtitle = "n = 500 | 50% censoring | dashed line: true break"
  )
  
  p_hist_h0 <- plot_split_points_hist(
    split_data[delta == 0], truth = tau, breaks = seq(-0.05, 10.05, by = 0.5),
    fill = pal_n[["500"]], x_breaks = 0:10, x_lab = expression(hat(tau)),
    subtitle = "False splits under H0 | n = 500 | 50% censoring"
  )
  
  ## 4. Save
  save_tables(list(
    invalid_jobs = invalid_jobs, censoring_check = censoring_check,
    summary_rejection = summary_rejection, summary_split_point = summary_split_point
  ), results_dir)
  
  save_plot(p_ecdf, "pvalues_ecdf_h0", plot_dir)
  save_plot(p_hist, "split_points_hist", plot_dir, width = 9, height = 6)
  save_plot(p_hist_h0, "split_points_hist_h0", plot_dir, width = 9, height = 4)
  EOF
  # test: write into scratchpad instead of results/
  sed -e "s#^results_dir <- results_path(run_id)#results_dir <- \"$S/pch_out\"#" $S/pch_analysis_skeleton.R > $S/pch_analysis_test.R
  echo 'print(censoring_check); print(summary_rejection); print(summary_split_point); print(invalid_jobs)' >> $S/pch_analysis_test.R
  cd /Users/jakobhaas/Documents/AAA_Bachelorarbeit/PAMOBs && Rscript $S/pch_analysis_test.R 2>&1 | grep -v "^Lade\|^Attach\|maskiert\|^$" | tail -25; ls $S/pch_out/figures
  
  OUT
  %notin%
    delta censoring_rate observed      diff
  <num>          <num>    <num>     <num>
    1:   0.0            0.5 0.500268  0.000268
  2:   0.5            0.5 0.499636 -0.000364
  3:   1.0            0.5 0.497970 -0.002030
  n censoring_rate delta n_valid n_splits split_rate   ci_lower  ci_upper
  <num>          <num> <num>   <int>    <int>      <num>      <num>     <num>
    1:   500            0.5   0.0    1000       58      0.058 0.04433286 0.0743365
  2:   500            0.5   0.5    1000      917      0.917 0.89814670 0.9333511
  3:   500            0.5   1.0    1000     1000      1.000 0.99631792 1.0000000
  n censoring_rate delta n_splits mean_split  sd_split       bias
  <num>          <num> <num>    <int>      <num>     <num>      <num>
    1:   500            0.5   0.5      917   5.023882 0.9461638 0.02388222
  2:   500            0.5   1.0     1000   5.011500 0.2015400 0.01150000
  mean_abs_error      rmse share_within_0.5
  <num>     <num>            <num>
    1:      0.5621592 0.9459493        0.6902944
  2:      0.1035000 0.2017672        0.9700000
  Empty data.table (0 rows and 3 cols): delta,n_invalid,share_invalid
  pvalues_ecdf_h0.pdf
  pvalues_ecdf_h0.png
  split_points_hist.pdf
  split_points_hist.png
  split_points_hist_h0.pdf
  split_points_hist_h0.png
  12:55 AM
  Das Skript läuft fehlerfrei durch. Ich prüfe noch kurz das Histogramm.
  
  
  
  Read split_points_hist.png
  12:55 AM
  Ich habe ein Analyse-Skript nach dem Muster deiner structural_change_hazard_main_analysis.R aufgebaut und getestet. Es läuft fehlerfrei durch, und die Grafiken sehen gut aus. Speicher es als simulation_batchtools/analysis/main/pemtree_piecewise_constant_hazard_analysis.R:
    
    
    ### Analysis: Piecewise constant hazard (split on time, correctly specified node model)
    source("simulation_batchtools/analysis/analysis_helpers.R")
  
  run_id <- "pemtree_piecewise_constant_hazard_v1"
  results_dir <- results_path(run_id)
  plot_dir <- file.path(results_dir, "figures")
  tau <- 5
  
  res <- load_results(run_id)
  mark_valid(res)
  
  ## 1. Checks
  invalid_jobs <- res[, .(n_invalid = sum(!valid), share_invalid = mean(!valid)), by = "delta"][n_invalid > 0]
  censoring_check <- check_censoring(res, by = c("delta", "censoring_rate"))
  
  ## 2. Type I error and power (one table, delta = 0 is the type I error)
  summary_rejection <- summarise_rejection(res[valid == TRUE], by = c("n", "censoring_rate", "delta"))
  
  p_ecdf <- plot_pvalue_ecdf(
    res[delta == 0 & valid == TRUE], group = "n",
    subtitle = "delta = 0 | n = 500 | 50% censoring"
  )
  
  ## 3. Split points (true tau = 5)
  split_data <- res[valid == TRUE & split == TRUE]
  summary_split_point <- summarise_split_points(
    split_data[delta > 0], truth = tau, by = c("n", "censoring_rate", "delta"), within = 0.5
  )
  
  p_hist <- plot_split_points_hist(
    split_data[delta > 0], truth = tau, breaks = seq(-0.05, 10.05, by = 0.1),
    facet = facet_thesis_wrap("delta", ncol = 1, scales = "free_y"),
    fill = pal_n[["500"]], x_breaks = 0:10, x_lab = expression(hat(tau)),
    subtitle = "n = 500 | 50% censoring | dashed line: true break"
  )

