#### setup
# packages
source("setup.R")
# functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemtree.R")
source("R/pemtree_splines.R")
# simulation helpers
source("simulation/R/pemtree_helpers.R")
source("simulation/R/run_one_level_shift_splines.R")
source("simulation/scenarios/02_level_shift_splines.R")

run_id <- "level_shift_null_splines_v1"
block_size <- 50

output_dir <- file.path("simulation/results/main", run_id, "blocks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

design <- expand.grid(
  n = c(1000, 2000),
  delta = 0,
  tau = 5,
  covariate_setting = c("none", "included"),
  df = 5,
  censoring_rate = c(0.2, 0.4, 0.6),
  rep = seq_len(1000),
  admin_time = 10,
  sim_interval = 0.05,
  ped_interval = 0.1,
  min_events = 25,
  alpha = 0.05,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design$seed <- 5555 + seq_len(nrow(design))

basehaz_fun <- function(t) -2.5 + 0.25 * t - 0.025 * t^2
covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -1, max = 1),
  x3 = list(distfun = rbinom, prob = 0.5, size = 1)
)
covariate_effects <- c(x1 = 0.3, x2 = 0.2, x3 = 0.5)


design$design_id <- seq_len(nrow(design))
design$block_id <- ceiling(design$design_id / block_size)

design_file <- file.path("simulation/results/main", run_id, "design.rds")
if (file.exists(design_file)) {
  saved_design <- readRDS(design_file)
  if (!identical(saved_design, design)) {
    stop("Current design differs from saved design")
  }
} else {
  saveRDS(design, design_file)
}

arg_names <- setdiff(colnames(design), c("design_id", "block_id"))

for (block in sort(unique(design$block_id))) {
  block_file <- file.path(output_dir, sprintf("block_%03d.rds", block))
  if (file.exists(block_file)) {
    message("Skipping existing block", block_file)
    next
  }

  rows <- which(design$block_id == block)
  message("Running block ", block, " of ", max(design$block_id))

  block_results <- do.call(
    rbind,
    lapply(rows, function(i) {
      args <- as.list(design[i, arg_names, drop = FALSE])
      args$basehaz_fun <- basehaz_fun
      if (args$covariate_setting == "included") {
        args$covariate_spec <- covariate_spec
        args$covariate_effects <- covariate_effects
      }
      do.call(run_one_level_shift_splines, args)
    })
  )

  block_results$design_id <- rows
  block_results$block_id <- block
  block_results$run_id <- run_id

  tmp_file <- paste0(block_file, ".tmp")
  saveRDS(block_results, tmp_file)

  if (!file.rename(tmp_file, block_file)) {
    stop("Could not finalise block ", block)
  }
}
