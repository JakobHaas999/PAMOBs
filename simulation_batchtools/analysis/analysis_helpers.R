# Shared helpers for all analysis scripts --------------------------------
# Source at the top of every analysis script (working directory = repo root):
#   source("simulation_batchtools/analysis/analysis_helpers.R")
#
# Everything that should look the same across the thesis (colours, shapes,
# labels, theme, axis formats) is defined here exactly once.

library("data.table")
library("ggplot2")


# Colours, shapes and labels used across all figures -------------------
pal_n <- c("250" = "#E69F00", "500" = "#009E73", "1000" = "#7B3294")
shape_n <- c("250" = 15, "500" = 16, "1000" = 17)
linetype_n <- c("250" = "solid", "500" = "longdash", "1000" = "dotdash")

linetype_df <- c("4" = "solid", "6" = "longdash", "8" = "dotdash")

pal_model <- c(oracle_log = "#000000", sp0 = "#0072B2", sp30 = "#D55E00", REML = "#CC79A7")
shape_model <- c(oracle_log = 18, sp0 = 15, sp30 = 16, REML = 17)
labels_model <- c(
  oracle_log = "Oracle: log(t)", sp0 = "P-spline, sp = 0",
  sp30 = "P-spline, sp = 30", REML = "P-spline, REML"
)

labels_scenario <- c(cubic = "Cubic baseline hazard", weibull = "Weibull baseline hazard")
labels_split_type <- c(binary = "Binary z", numeric = "Numeric z")
labels_covariate_setting <- c(none = "Without covariates", included = "With covariates")
label_censoring <- function(x) paste0(100 * as.numeric(x), "% censoring")
label_df <- function(x) paste0("df = ", x)
label_delta <- function(x) paste0("delta = ", x)
label_percent_x <- function(x) scales::percent(as.numeric(x), accuracy = 1)
labels_tree_model <- c(pemtree = "PEM tree", weibull_aft = "Weibull AFT tree")

# One labeller for all facets; unknown variables keep their raw value
thesis_labeller <- labeller(
  scenario = labels_scenario,
  split_type = labels_split_type,
  model = labels_tree_model,
  covariate_setting = labels_covariate_setting,
  censoring_rate = label_censoring,
  df = label_df,
  delta = label_delta,
  .default = label_value
)


# Scales, facets and theme ------------------------------------------------
# Ready-made scale bundles (a list can be added to a ggplot with +)
scale_n <- function(aesthetics = c("colour", "shape")) {
  list(
    colour = scale_colour_manual(values = pal_n, name = "Sample size"),
    shape = scale_shape_manual(values = shape_n, name = "Sample size"),
    fill = scale_fill_manual(values = pal_n, name = "Sample size")
  )[aesthetics]
}

scale_model <- function() {
  list(
    scale_colour_manual(values = pal_model, labels = labels_model, name = "Node model"),
    scale_shape_manual(values = shape_model, labels = labels_model, name = "Node model")
  )
}

# facet_grid / facet_wrap with the thesis labeller
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
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
  if (discrete_x) th <- th + theme(panel.grid.major.x = element_blank())
  th
}


# Paths and loading -------------------------------------------------------
results_path <- function(run_id, phase = "main") {
  file.path("simulation_batchtools/results", phase, run_id)
}

registry_path <- function(run_id, phase = "main") {
  file.path("simulation_batchtools/registries", phase, run_id)
}

# Load job results: from the cached job_results.rds if available, otherwise
# from the batchtools registry (and write the cache). Jobs that return a
# table with several rows (one per node model) are unpacked automatically.
load_results <- function(run_id, phase = "main") {
  cache <- file.path(results_path(run_id, phase), "job_results.rds")
  if (file.exists(cache)) {
    return(as.data.table(readRDS(cache)))
  }

  reg <- suppressWarnings(batchtools::loadRegistry(
    registry_path(run_id, phase),
    work.dir = getwd(), writeable = FALSE, make.default = FALSE
  ))
  res <- batchtools::reduceResultsDataTable(reg = reg)
  if (is.data.frame(res$result[[1]])) {
    res <- rbindlist(
      Map(function(id, x) cbind(job.id = id, as.data.table(x)), res$job.id, res$result),
      fill = TRUE
    )
  } else {
    res <- batchtools::unwrap(res)
  }
  results <- merge(batchtools::unwrap(batchtools::getJobPars(reg = reg)), res, by = "job.id")

  dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
  saveRDS(results, cache)
  results
}

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
plot_rejection_by_censoring <- function(data, group = c("n", "variant"), facet = NULL,
                                        y_max = NULL, ref_line = 0.05,
                                        title = NULL, subtitle = NULL,
                                        y_lab = "Rejection Rate") {
  group <- match.arg(group)
  levels_group <- if (group == "n") names(pal_n) else names(pal_model)
  dodge <- position_dodge(width = 0.35)

  y_scale <- if (is.null(y_max)) {
    scale_y_continuous(
      limits = c(0, NA), breaks = scales::breaks_pretty(),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.1))
    )
  } else {
    scale_y_continuous(
      limits = c(0, y_max), breaks = seq(0, y_max, by = 0.2),
      labels = scales::label_percent(accuracy = 1),
      expand = expansion(mult = c(0.02, 0.03))
    )
  }

  ggplot(
    data,
    aes(
      x = factor(censoring_rate), y = split_rate,
      colour = factor(.data[[group]], levels = levels_group),
      shape = factor(.data[[group]], levels = levels_group),
      group = factor(.data[[group]], levels = levels_group)
    )
  ) +
    (if (!is.null(ref_line)) {
      geom_hline(yintercept = ref_line, linetype = "dashed", colour = "grey45", linewidth = 0.5)
    }) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = dodge, width = 0.1, linewidth = 0.6) +
    geom_point(position = dodge, size = 3) +
    facet +
    (if (group == "n") scale_n() else scale_model()) +
    scale_x_discrete(labels = label_percent_x) +
    y_scale +
    labs(title = title, subtitle = subtitle, x = "Target censoring rate", y = y_lab) +
    theme_thesis(discrete_x = TRUE)
}

# Power curve: split rate against the effect size, one line per sample size.
plot_power_curve <- function(data, x = "delta", facet = NULL, x_lab = expression(delta),
                             x_breaks = waiver(), title = NULL, subtitle = NULL,
                             y_lab = "Rejection Rate") {
  ebar_width <- 0.03 * diff(range(data[[x]]))
  ggplot(
    data,
    aes(
      x = .data[[x]], y = split_rate,
      colour = factor(n, levels = names(pal_n)),
      shape = factor(n, levels = names(pal_n)),
      group = factor(n, levels = names(pal_n))
    )
  ) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
    geom_line(linewidth = 0.6) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = ebar_width, linewidth = 0.5) +
    geom_point(size = 2.5) +
    facet +
    scale_n() +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(
      limits = c(0, 1), breaks = seq(0, 1, by = 0.2),
      labels = scales::label_percent(accuracy = 1)
    ) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab) +
    theme_thesis()
}

# Boxplots of estimated split points, filled by sample size.
plot_split_points_box <- function(data, truth, x = "delta", split_col = "split_point",
                                  facet = NULL, y_limits = NULL, y_breaks = waiver(),
                                  x_lab = expression(delta), x_labels = waiver(),
                                  y_lab = "Estimated split point",
                                  title = NULL, subtitle = NULL) {
  ggplot(
    data,
    aes(
      x = factor(.data[[x]]), y = .data[[split_col]],
      fill = factor(n, levels = names(pal_n))
    )
  ) +
    geom_hline(yintercept = truth, linetype = "dashed", colour = "grey45", linewidth = 0.5) +
    geom_boxplot(outlier.size = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
    facet +
    scale_n("fill") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = y_lab) +
    theme_thesis(discrete_x = TRUE)
}

# Histogram of estimated split points; pass bin breaks centred on the grid.
plot_split_points_hist <- function(data, truth, breaks, split_col = "split_point",
                                   facet = NULL, x_limits = NULL, x_breaks = waiver(),
                                   fill = pal_n[["1000"]], x_lab = "Estimated split point",
                                   title = NULL, subtitle = NULL) {
  ggplot(data, aes(x = .data[[split_col]])) +
    geom_histogram(breaks = breaks, fill = fill, colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = truth, linetype = "dashed", colour = "grey45") +
    facet +
    scale_x_continuous(breaks = x_breaks) +
    coord_cartesian(xlim = x_limits) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = "Count") +
    theme_thesis()
}

# ECDF of root-test p-values (should follow the diagonal under H0).
#   group = "n"  -> colour + linetype by sample size
#   group = "df" -> linetype by spline df
plot_pvalue_ecdf <- function(data, group = c("n", "df"), p_col = "p_root", facet = NULL,
                             title = NULL, subtitle = NULL) {
  group <- match.arg(group)
  mapping <- if (group == "n") {
    aes(
      x = .data[[p_col]],
      colour = factor(n, levels = names(pal_n)),
      linetype = factor(n, levels = names(pal_n))
    )
  } else {
    aes(x = .data[[p_col]], linetype = factor(df, levels = names(linetype_df)))
  }
  group_scales <- if (group == "n") {
    list(
      scale_colour_manual(values = pal_n, name = "Sample size"),
      scale_linetype_manual(values = linetype_n, name = "Sample size")
    )
  } else {
    scale_linetype_manual(values = linetype_df, name = "Spline df")
  }

  ggplot(data, mapping) +
    geom_abline(intercept = 0, slope = 1, colour = "grey60", linetype = "dashed") +
    stat_ecdf(geom = "step", linewidth = 0.8) +
    facet +
    group_scales +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
    labs(title = title, subtitle = subtitle, x = "Root p-value", y = "Cumulative proportion") +
    theme_thesis()
}


# Saving ------------------------------------------------------------------
# Save a plot as PDF (for LaTeX) and high-resolution PNG
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
