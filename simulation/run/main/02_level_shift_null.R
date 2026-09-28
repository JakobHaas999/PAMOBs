### source relevant files
## packages
source("setup.R")

## functions
source("R/sim_pexp_own.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemtree.R")

## simulation helpers
source("simulation/R/fit_models.R")
source("simulation/R/run_one.R")

## simulation scenario
source("simulation/scenarios/01_level_shift.R")

run_id <- "level_shift_null_v1"
block_size <- 50

output_dir <- file.path("simulation/results/main", run_id, "blocks")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

design <- expand.grid(
  n = c(500, 1000),
  delta = 0,
  tau = 4,
  beta0 = -2,
  censoring_rate = c(0.3, 0.5, 0.7),
  rep = seq_len(500),
  ped_interval = 0.1,
  alpha = 0.05,
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 200000 + seq_len(nrow(design))

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
    message("Skipping existing block", block)
    next
  }

  rows <- which(design$block_id == block)
  message("Running block ", block, " of ", max(design$block_id))

  block_results <- do.call(
    rbind,
    lapply(rows, function(i) {
      do.call(
        run_one_level_shift,
        as.list(design[i, arg_names, drop = FALSE])
      )
    })
  )

  block_results$design_id <- rows
  block_results$block_id <- block
  block_results$run_id <- run_id

  tmp_file <- paste0(block_file, ".tmp")
  saveRDS(block_results, tmp_file)

  if (!file.rename(tmp_file, block_file)) {
    stop("Could not finalize block ", block)
  }
}
