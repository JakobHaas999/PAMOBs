# -------------------------------------------------------------
# Function to read simulation result block data
# -------------------------------------------------------------

read_block_data <- function(run_id) {
  run_dir <- file.path("simulation/results/main", run_id)
  blocks_dir <- file.path(run_dir, "blocks")

  block_files <- list.files(blocks_dir, pattern = "^block_[0-9]+\\.rds$", full.names = TRUE)

  if (!length(block_files)) stop("No result blocks found for run_id: ", run_id)

  data.table::rbindlist(
    lapply(block_files, readRDS),
    use.names = TRUE, fill = TRUE
  )
}

# -------------------------------------------------------------
# Function to summarize empirical type I error rates
# -------------------------------------------------------------

summary_type_I_error <- function(results, grouping_cols = c("scenario", "n", "censoring_rate")) {
  summary <- results[,
    {
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
          ci_low = test$conf.int[[1]],
          ci_high = test$conf.int[[2]],
          p_value = test$p.value
        )
      }
    },
    by = grouping_cols
  ]

  data.table::setorder(summary, scenario, n, censoring_rate)
  return(summary)
}

# -------------------------------------------------------------
# Function to summarize observed censoring rates
# -------------------------------------------------------------

summary_obs_censoring <- function(results, grouping_cols) {
  summary <- results[success %in% TRUE, .(
    mean_observed_censoring = mean(observed_censoring),
    sd_observed_censoring = sd(observed_censoring)
  ), by = grouping_cols]

  data.table::setorder(summary, scenario, n, censoring_rate)
  return(summary)
}


# -------------------------------------------------------------
# Function to plot type-I-error estimates
# -------------------------------------------------------------

plot_type_I_error <- function(summary) {
  plot_data <- data.table::copy(summary)
  plot_data[, n_label := factor(n)]

  scenario_labels <- c(
    "level shift with covariates" = "With covariates",
    "level shift without covariates" = "Without covariates"
  )

  point_position <- position_dodge(width = .025)
  y_upper <- max(.06, 1.1 * max(plot_data$ci_high, na.rm = TRUE))

  type1_plot <- ggplot(plot_data, aes(
    x = censoring_rate, y = alpha_hat,
    colour = n_label,
    shape = n_label,
    group = n_label
  )) +
    geom_hline(yintercept = .05, linetype = "dashed", colour = "grey35", linewidth = .7) +
    geom_line(position = point_position, linewidth = .7) +
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high),
      width = .015, linewidth = .65
    ) +
    geom_point(position = point_position, size = 2.8) +
    facet_wrap(~scenario, labeller = as_labeller(scenario_labels)) +
    scale_x_continuous(breaks = c(.2, .4, .6), labels = scales::label_percent(accuracy = 1)) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 0.1), expand = expansion(mult = c(0, .04))) +
    scale_colour_manual(values = c("1000" = "#0072B2", "2000" = "#D55E00")) +
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
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold")
    )

  print(type1_plot)
}
