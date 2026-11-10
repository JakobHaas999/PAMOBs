### Analysis of complex tree (two or more splits on z and time). -----------------
source("simulation_batchtools/analysis/analysis_helpers.R")

run_id <- "complex_tree_v1"
results_dir <- results_path(run_id)
plot_dir <- file.path(results_dir, "figures")

results <- load_results(run_id)
mark_valid(results)

results[, structure := fcase(
  root_var == "z" & right_var == "tend" & left_var == "none" & n_split == 2, "Correct tree",
  root_var == "z" & right_var == "tend", "Correct tree + extra splits",
  root_var == "z" & right_var == "none" & left_var == "none", "Time split missed",
  default = "Other"
)]
results[, structure := factor(structure, levels = c(
  "Correct tree", "Correct tree + extra splits", "Time split missed", "Other"
))]

share <- results[, .N, by = c("scenario", "dz", "dt", "censoring_rate", "structure")]
share[, prop := N / sum(N), by = c("scenario", "dz", "dt", "censoring_rate")]

p_structure_baseline_shift <- ggplot(
  share[scenario == "covar_struc_chang"],
  aes(x = factor(censoring_rate), y = prop, fill = structure)
) +
  geom_col(width = 0.7) +
  facet_thesis(rows = "dt", cols = "dz") +
  scale_fill_manual(values = pal_structure) +
  scale_x_discrete(labels = label_percent_x) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "Recovered tree structure | Baseline shift in the older group",
    subtitle = "True tree: split on z at 0, then on time at 3 in the older branch only | n = 1000",
    x = "Target censoring rate", y = "Share of trees", fill = NULL
  ) +
  theme_thesis(discrete_x = TRUE)
p_structure_baseline_shift

p_structure_time_chang <- ggplot(
  share[scenario == "covar_time_chang"],
  aes(
    x = factor(censoring_rate),
    y = prop,
    fill = structure
  )
) +
  geom_col(width = 0.7) +
  facet_thesis(rows = "dt", cols = "dz") +
  scale_fill_manual(values = pal_structure) +
  scale_x_discrete(labels = label_percent_x) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    title = "Recovered tree structure | Time-varying treatment effect with z > 0",
    subtitle = "True tree: split on z at 0, then on time at 3 in the group with z > 0 only | n = 1000",
    x = "Target censoring rate", y = "Share of trees", fill = NULL
  ) +
  theme_thesis(discrete_x = TRUE)
p_structure_time_chang
