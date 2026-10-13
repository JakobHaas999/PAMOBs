########################################################
## Analysis for 03 Comparison Pemtree Pemob  ###########
########################################################

### setup
source("setup.R")
source("simulation/R/analysis_helpers.R")

run_id <- "compare_pemob_pemtree_v1"
results <- read_block_data(run_id)

required_columns <- c(
  "design_id", "n", "censoring_rate", "success",
  "split_detected", "observed_censoring", "method"
)
missing_columns <- setdiff(required_columns, colnames(results))
if (length(missing_columns)) stop("Missing result columns", paste(missing_columns, ", "))

data.table::setorder(results, design_id)

### Type-I error summary
size_summary <- summary_type_I_error(results, grouping_cols = c("method", "n", "censoring_rate"))
data.table::setorder(size_summary, method)
size_summary

### Observed censoring summary
censoring_summary <- summary_obs_censoring(results, grouping_cols = c("method", "n", "censoring_rate"))
data.table::setorder(censoring_summary, method)
censoring_summary

### Plot summary
plot_data <- copy(size_summary)[, n_label := factor(n)]
point_position <- position_dodge(width = .025)
y_upper <- max(.06, 1.1 * max(plot_data$ci_high, na.rm = TRUE))

type1_plot <- ggplot(plot_data, aes(
  x = censoring_rate, y = alpha_hat,
  colour = n_label, group = n_label, shape = n_label
)) +
  geom_hline(yintercept = .05, linetype = "dashed", colour = "grey35", linewidth = .7) +
  geom_line(position = point_position, linewidth = .7) +
  geom_point(position = point_position, size = 2.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high)) +
  facet_wrap(~method) +
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
