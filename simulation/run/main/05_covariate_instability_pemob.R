### setup
source("setup.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemob.R")
# Simulation helper
source("simulation/R/run_covariate_instability.R")

### Run information
run_id <- "covariate_instability_pemob_v1"
block_size <- 50

### Specify simulation parameters

# Quadratic baseline hazard
log_hazard <- function(t) {
  -3 + 0.15 * t - 0.01 * t^2
}

# Covariate specifications
covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = 0, max = 3),
  x3 = list(distfun = rbinom, prob = 0.5, size = 1)
)

# Formula for simulation
sim_formula <- function(d) {
  reformulate(
    termlabels = c(
      "log_hazard(t)",
      "0.7 * x1",
      "0.4 * x2",
      sprintf("(%s) * (x3 == 1) * x1", d$delta)
    ), env = environment()
  )
}
fit_formula <- Surv(time, status) ~ bs(tend, df = 5) + x1 + x2 + offset(offset) | x3

# Design with parameters
design <- expand.grid(
  n = c(1000, 2000),
  delta = seq(0, 1, by = 0.2),
  censoring_rate = c(0.2, 0.4, 0.6),
  rep = seq_len(1000),
  admin_time = 10,
  sim_interval = 0.05,
  fit_interval = 0.1,
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 10000 + seq_len(nrow(design))

design$design_id <- seq_len(nrow(design))
design$block_id <- ceiling(design$design_id / block_size)

### Output
output_dir <- file.path("simulation/results/main", run_id)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

### Save design
design_file <- file.path(output_dir, "design.rds")
if (file.exists(design_file)) {
  saved_design <- readRDS(design_file)
  if (!identical(saved_design, design)) {
    stop("Current design differs from saved design")
  }
} else {
  saveRDS(design, design_file)
}

### Run simulations per block
for (block in unique(sort(design$block_id))) {
  block_file <- file.path(output_dir, sprintf("block_%03d.rds", block))
  if (file.exists(block_file)) {
    message("Skipping existing block: ", block)
    next
  }
  rows <- which(design$block_id == block)
  message("Running block ", block, " of ", max(design$block_id))
  block_design <- design[rows, , drop = FALSE]
  block_results <- run_covariate_instability(
    log_hazard = log_hazard,
    design = block_design,
    sim_formula = sim_formula,
    fit_formula = fit_formula,
    fit_fun = pemob,
    covariate_spec = covariate_spec,
    alpha = 0.05,
    maxdepth = 2
  )
  block_results$run_id <- run_id
  tmp_file <- paste0(block_file, ".tmp")
  saveRDS(block_results, tmp_file)
  if (!file.rename(tmp_file, block_file)) {
    stop("Could not finalise block", block)
  }
}
