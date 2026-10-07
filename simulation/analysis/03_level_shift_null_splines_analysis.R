####################################################################
## Analysis for 03 level shift null scenario with splines ##########
####################################################################

### setup
source("setup.R")

run_id <- "level_shift_null_splines_v1"
output_dir <- file.path("simulation/results/main", run_id, "blocks")
files <- list.files(output_dir)

results <- rbindlist(lapply(files, function(f) {
  as.data.table(readRDS(file.path(output_dir, f)))
}))

## create results summary with alpha_hat, mcse and CI per scenario, n, censoring_rate
size_summary <- results[,
  {
    valid <- success %in% TRUE & !is.na(split_detected)

    n_success <- sum(valid)
    n_splits <- sum(split_detected[valid])

    if (n_success == 0L) {
      list(
        n_runs = .N, n_success = 0, n_errors = .N, n_splits = NA_integer_,
        alpha_hat = NA_real_, mcse = NA_real_, ci_low = NA_real_, ci_high = NA_real_
      )
    } else {
      alpha_hat <- n_splits / n_success
      test <- binom.test(
        x = n_splits,
        n = n_success,
        p = 0.05,
        conf.level = 0.95
      )

      list(
        n_runs = .N, n_success = n_success, n_errors = .N - n_success, n_splits = n_splits,
        alpha_hat = alpha_hat, mcse = sqrt(alpha_hat * (1 - alpha_hat) / n_success),
        ci_low = unname(test$conf.int[1]), ci_high = unname(test$conf.int[2])
      )
    }
  },
  by = .(scenario, n, censoring_rate)
]
setorder(size_summary, scenario, n, censoring_rate)
size_summary

## Plot results
plot_data <- data.table::copy(size_summary)
plot_data[, n_label := factor(n)]

scenario_labels <- c(
  "level shift with covariates" = "With covariates",
  "level shift without covariates" = "Without covariates"
)

type1_plot <- ggplot(
  plot_data,
  aes(x = censoring_rate, y = alpha_hat, colour = n_label, group = n_label)
) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey80", linewidth = .7) +
  geom_line(linewidth = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .015, linewidth = .65) +
  geom_point(size = 2.8) +
  facet_wrap(~scenario, labeller = as_labeller(scenario_labels)) +
  scale_x_continuous(breaks = c(.2, .4, .6), labels = scales::label_percent(accuracy = 1)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = .1), expand = expansion(mult = c(0, .8))) +
  scale_colour_manual(values = c("1000" = "#0072B2", "2000" = "#D55E00")) +
  labs(
    x = "Censoring rate",
    y = "Empirical Type-I-Error",
    colour = "sample size"
  ) +
  theme_bw() +
  theme(
    legend.position = "top"
  )
type1_plot

ggsave("simulation/analysis/figures/type1_error_splines.png",
  plot = type1_plot,
  width = 8,
  height = 5,
  dpi = 300
)
