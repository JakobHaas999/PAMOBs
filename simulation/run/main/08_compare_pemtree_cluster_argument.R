### setup
# Packages
source("setup.R")
# Functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemtree_splines.R")

##### Main simulation

### Run info
run_id <- "compare_pemtree_cluster_argument_v1"
block_size <- 50

### Parameters
log_hazard <- function(t, delta) {
  -2 + delta * (t >= 5)
}

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -1, max = 1)
)

design <- expand.grid(
  n = c(1000, 2000),
  delta = c(0, 0.5),
  censoring_rate = c(0.2, 0.4, 0.6),
  covariate_setting = c("none", "include"),
  rep = seq_len(1000),
  admin_time = 10,
  sim_interval = 0.05,
  ped_interval = 0.1,
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 600000 + seq_len(nrow(design))
design$design_id <- seq_len(nrow(design))
design$block_id <- ceiling(design$design_id / block_size)

### Output
output_dir <- file.path("simulation/results/main", run_id)
blocks_dir <- file.path(output_dir, "blocks")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(blocks_dir, showWarnings = FALSE, recursive = TRUE)

### Save design
design_file <- file.path(output_dir, "design.rds")
if (file.exists(design_file)) {
  saved_design <- readRDS(design_file)
  if (!identical(saved_design, design)) {
    stop("Saved design differs from saved design")
  }
} else {
  saveRDS(design, design_file)
}

### Run simulatiom
for (block in unique(sort(design$block_id))) {
  block_file <- file.path(blocks_dir, sprintf("block_%03d.rds", block))
  if (file.exists(block_file)) {
    message("Skipping block: ", block)
    next
  }
  rows <- which(design$block_id == block)
  message("Running block ", block, " of ", max(design$block_id))

  block_design <- design[rows, , drop = FALSE]

  block_results <- lapply(seq_len(nrow(block_design)), function(i) {
    d <- block_design[i, , drop = FALSE]
    set.seed(d$seed)

    terms <- sprintf("log_hazard(t, delta = %s)", d$delta)
    if (d$covariate_setting == "none") {
      fit_formula <- Surv(time, status) ~ offset(offset) | tend
      covar_spec <- NULL
    } else {
      fit_formula <- Surv(time, status) ~ offset(offset) + x1 + x2 | tend
      covar_spec <- covariate_spec
      terms <- c(terms, "0.7 * x1", "0.4 * x2")
    }
    sim_formula <- reformulate(terms)

    out <- tryCatch(
      {
        sim_data <- sim_pam(
          n = d$n,
          formula = sim_formula,
          covariate_spec = covar_spec,
          censoring_rate = d$censoring_rate,
          admin_time = d$admin_time,
          cut = seq(0, d$admin_time, by = d$sim_interval)
        )

        args <- list(
          formula = fit_formula,
          data = sim_data,
          cut = seq(0, d$admin_time, by = d$ped_interval),
          alpha = 0.05,
          maxdepth = 2
        )

        fit_normal <- do.call(pemtree_splines, c(args, list(cluster_id = NULL)))
        fit_cluster <- do.call(pemtree_splines, c(args, list(cluster_id = "id")))

        fits <- list(normal = fit_normal, cluster = fit_cluster)
        tree_out <- Map(
          function(fit, method) {
            n_terminal <- length(nodeids(fit$tree, terminal = TRUE))
            split_detected <- n_terminal > 1

            tau_hat <- if (split_detected) {
              root <- node_party(fit$tree)
              sp <- split_node(root)
              sp$breaks[1]
            } else {
              NA_real_
            }

            cbind(
              d,
              method = method,
              success = TRUE,
              n_terminal_nodes = n_terminal,
              split_detected = split_detected,
              tau_hat = tau_hat,
              observed_censoring = mean(sim_data$status == 0),
              error = NA_character_
            )
          }, fits, names(fits)
        )
        do.call(rbind, tree_out)
      },
      error = function(e) {
        cbind(
          d,
          method = c("normal", "cluster"),
          success = FALSE,
          n_terminal_nodes = NA_integer_,
          split_detected = NA,
          tau_hat = NA_real_,
          observed_censoring = NA_real_,
          error = conditionMessage(e)
        )
      }
    )
    out
  })
  block_results <- do.call(rbind, block_results)
  rownames(block_results) <- NULL

  saveRDS(block_results, block_file)
  message("Saved block ", block)
}
