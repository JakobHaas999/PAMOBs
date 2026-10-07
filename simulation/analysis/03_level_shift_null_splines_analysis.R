####################################################################
## Analysis for 03 level shift null scenario with splines ##########
####################################################################

### setup
source("setup.R")

run_id <- "level_shift_null_splines_v1"
run_dir <- file.path("simulation/results/main", run_id)
blocks_dir <- file.path(run_dir, "blocks")
design_file <- file.path(run_dir, "design.rds")

block_files <- list.files(
  blocks_dir,
  pattern = "^block_[0-9]+\\.rds$",
  full.names = TRUE
)

if (length(block_files) == 0L) {
  stop("No result blocks found for run_id: ", run_id)
}

results <- data.table::rbindlist(
  lapply(block_files, readRDS),
  use.names = TRUE,
  fill = TRUE
)

required_columns <- c(
  "design_id", "scenario", "n", "censoring_rate", "success",
  "split_detected", "observed_censoring"
)
missing_columns <- setdiff(required_columns, names(results))
if (length(missing_columns) > 0L) {
  stop("Missing result columns: ", paste(missing_columns, collapse = ", "))
}

if (anyDuplicated(results$design_id)) {
  stop("Duplicate design_id values found in the result blocks")
}

if (file.exists(design_file)) {
  design <- readRDS(design_file)
  missing_ids <- setdiff(design$design_id, results$design_id)
  unexpected_ids <- setdiff(results$design_id, design$design_id)

  if (length(missing_ids) > 0L || length(unexpected_ids) > 0L) {
    stop(
      "Results do not match the saved design: ",
      length(missing_ids), " missing and ",
      length(unexpected_ids), " unexpected design IDs"
    )
  }
}

data.table::setorder(results, design_id)

### Type-I error summary
size_summary <- results[, {
  valid <- success %in% TRUE & !is.na(split_detected)
  n_success <- sum(valid)
  n_splits <- sum(split_detected[valid])

  if (n_success == 0L) {
    list(
      n_runs = .N,
      n_success = 0L,
      n_errors = .N,
      n_splits = NA_integer_,
      alpha_hat = NA_real_,
      mcse = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_
    )
  } else {
    alpha_hat <- n_splits / n_success
    test <- stats::binom.test(
      x = n_splits,
      n = n_success,
      p = 0.05,
      conf.level = 0.95
    )

    list(
      n_runs = .N,
      n_success = n_success,
      n_errors = .N - n_success,
      n_splits = n_splits,
      alpha_hat = alpha_hat,
      mcse = sqrt(alpha_hat * (1 - alpha_hat) / n_success),
      ci_low = unname(test$conf.int[1]),
      ci_high = unname(test$conf.int[2]),
      p_value = test$p.value
    )
  }
}, by = .(scenario, n, censoring_rate)]

data.table::setorder(size_summary, scenario, n, censoring_rate)
print(size_summary)

### Observed censoring summary
censoring_summary <- results[success %in% TRUE, .(
  mean_observed_censoring = mean(observed_censoring),
  sd_observed_censoring = stats::sd(observed_censoring)
), by = .(scenario, n, censoring_rate)]

data.table::setorder(censoring_summary, scenario, n, censoring_rate)
print(censoring_summary)

### Plot Type-I error estimates
plot_data <- data.table::copy(size_summary)
plot_data[, n_label := factor(n)]

scenario_labels <- c(
  "level shift with covariates" = "With covariates",
  "level shift without covariates" = "Without covariates"
)

point_position <- position_dodge(width = 0.025)
y_upper <- max(0.06, 1.1 * max(plot_data$ci_high, na.rm = TRUE))

type1_plot <- ggplot(
  plot_data,
  aes(
    x = censoring_rate,
    y = alpha_hat,
    colour = n_label,
    shape = n_label,
    group = n_label
  )
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    colour = "grey35",
    linewidth = 0.7
  ) +
  geom_line(position = point_position, linewidth = 0.7) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = point_position,
    width = 0.015,
    linewidth = 0.65
  ) +
  geom_point(position = point_position, size = 2.8) +
  facet_wrap(
    ~scenario,
    labeller = as_labeller(scenario_labels)
  ) +
  scale_x_continuous(
    breaks = c(0.2, 0.4, 0.6),
    labels = scales::label_percent(accuracy = 1)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.04))
  ) +
  scale_colour_manual(
    values = c("1000" = "#0072B2", "2000" = "#D55E00")
  ) +
  scale_shape_manual(values = c("1000" = 16, "2000" = 17)) +
  coord_cartesian(ylim = c(0, y_upper)) +
  labs(
    x = "Censoring rate",
    y = "Empirical Type-I error",
    colour = "Sample size",
    shape = "Sample size"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold")
  )

print(type1_plot)

figures_dir <- file.path("simulation/analysis", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(figures_dir, "type1_error_splines.png"),
  plot = type1_plot,
  width = 8,
  height = 5,
  dpi = 300
)
