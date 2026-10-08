########################################################
## Analysis for 03 level shift null scenario  ##########
########################################################

### setup
source("setup.R")
source("simulation/R/analysis_helpers.R")

run_id <- "level_shift_null_constant_v1"
results <- read_block_data(run_id)

required_columns <- c(
  "design_id", "scenario", "n", "censoring_rate", "success",
  "split_detected", "observed_censoring"
)
missing_columns <- setdiff(required_columns, colnames(results))
if (length(missing_columns)) stop("Missing result columns", paste(missing_columns, ", "))
if (anyDuplicated(results$design_id)) stop("Duplicate design_id values found in the result blocks")

data.table::setorder(results, design_id)

### Type-I error summary
size_summary <- summary_type_I_error(results, c("scenario", "n", "censoring_rate"))
print(size_summary)

### Observed censoring summary
censoring_summary <- summary_obs_censoring(results, c("scenario", "n", "censoring_rate"))
print(censoring_summary)

### Plot Type-I error estimates
type1_plot <- plot_type_I_error(size_summary)

figures_dir <- file.path("simulation/analysis", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = file.path(figures_dir, "type1_error_constant.png"),
  plot = type1_plot,
  width = 8,
  height = 5,
  dpi = 300
)
