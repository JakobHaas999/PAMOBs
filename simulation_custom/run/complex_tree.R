##########################################################################
## Complex tree DGP: can pemtree reconstruct a tree with 3 split levels? ##

### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
library("parallel")

n_sim <- 10
n_obs <- 1000
n_pilot <- 5000

admin_time <- 10
censoring_rate <- 0.4

alpha <- 0.05
maxdepth <- 5
cut <- seq(0, admin_time, by = 0.1)

n_cores <- 6

## Baseline
shape <- 1.5
rate <- 0.05
log_baseline_hazard <- function(t) {
  log(shape) + log(rate) + (shape - 1) * log(t)
}

# Leaf assignments and leaf specified effects
true_leaf <- function(z1, z2, z3) {
  ifelse(
    z1 <= 0,
    ifelse(z2 <= 6, 1, ifelse(z3 <= 10, 2, 3)),
    ifelse(z3 <= 10, 4, 5)
  )
}
b1 <- c(0.5, 1, 1, -0.5, -0.5)
b2 <- c(0.5, -0.5, 0.5, 0, 1)

sim_formula <- ~ log_baseline_hazard(t) +
  b1[true_leaf(z1, z2, z3)] * x1 +
  b2[true_leaf(z1, z2, z3)] * x2

true_splits <- data.frame(
  path = c("", "L", "LR", "R"),
  variable = c("z1", "z2", "z3", "z3"),
  break_point = c(0, 6, 10, 10)
)

# Function to generate data
gen_data <- function(n) {
  d <- data.frame(
    id = seq_len(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    z1 = rnorm(n),
    z2 = rnorm(n, mean = 6, sd = 3),
    z3 = rnorm(n, mean = 10, sd = 2),
    z4 = rnorm(n, mean = 0, sd = 2),
    z5 = rnorm(n, mean = 0, sd = 2),
    z6 = runif(n, min = -1, max = 2)
  )
  d$true_leaf <- true_leaf(d$z1, d$z2, d$z3)
  d
}

fit_formula <- Surv(time, status) ~ bs(tend, df = 6) + x1 + x2 + offset(offset) |
  z1 + z2 + z3 + z4 + z5 + z6

extract_splits <- function(node, var_names, path = "") {
  if (is.terminal(node)) {
    return(NULL)
  }
  split <- split_node(node)
  kids <- kids_node(node)
  rbind(
    data.frame(
      path = path,
      variable = var_names[varid_split(split)],
      break_point = breaks_split(split)
    ),
    extract_splits(kids[[1]], var_names, paste0(path, "L")),
    extract_splits(kids[[2]], var_names, paste0(path, "R"))
  )
}

# Calibrate censoring
set.seed(123)
pilot_events <- sim_event_times(
  formula = sim_formula,
  data = gen_data(n_pilot),
  maxt = admin_time
)
lambda_c <- calibrate_censoring(
  event_times = pilot_events,
  censoring_rate = censoring_rate,
  admin_time = admin_time
)

# Simulation
seeds <- 7100 + seq_len(n_sim)

run_rep <- function(i) {
  message("Simulation ", i, " of", n_sim, " started")
  set.seed(seeds[i])

  sim_data <- sim_surv(
    formula = sim_formula,
    data = gen_data(n_obs),
    admin_time = admin_time,
    lambda_c = lambda_c
  )

  fit <- tryCatch(
    pemtree(
      formula = fit_formula,
      data = sim_data,
      cut = cut,
      maxdepth = maxdepth,
      alpha = alpha
    ),
    error = function(e) e
  )

  success <- !inherits(fit, "error")

  summary <- data.frame(
    rep = i,
    seed = seeds[[i]],
    success = success,
    error = if (success) NA_character_ else conditionMessage(fit),
    n_leaves = NA_integer_,
    exact_structure = NA,
    noise_selected = NA,
    observed_censoring = mean(sim_data$status == 0)
  )
  splits <- NULL
  nodes <- NULL

  if (success) {
    tree <- fit$tree
    splits <- extract_splits(node_party(tree), names(tree$data))
    if (is.null(splits)) {
      splits <- data.frame(path = character(0), variable = character(0), break_point = character(0))
    }

    nodes <- merge(true_splits, splits, by = "path", suffixes = c("_true", "_est"), all.x = TRUE)
    nodes$found <- !is.na(nodes$variable_est)
    nodes$variable_correct <- nodes$found & nodes$variable_est == nodes$variable_true
    nodes$break_error <- ifelse(nodes$variable_correct, nodes$break_point_est - nodes$break_point_true, NA_real_)

    summary$n_leaves <- length(nodeids(tree, terminal = TRUE))
    summary$exact_structure <- nrow(splits) == nrow(true_splits) && all(nodes$variable_correct)
    summary$noise_selected <- any(splits$variable %in% c("z4", "z5", "z6"))

    splits$rep <- i
    nodes$rep <- i
  }

  file.create(file.path(progress_dir, i))
  n_done <- length(list.files(progress_dir))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf(
    "[%s] Simulation %d done | %d/%d (%.0f%%) | %.1f min elapsed, ~%.1f min left",
    format(Sys.time(), "%H:%M:%S"), i, n_done, n_sim, 100 * n_done / n_sim,
    elapsed, elapsed / n_done * (n_sim - n_done)
  ))

  list(summary = summary, splits = splits, nodes = nodes)
}

progress_dir <- file.path(tempdir(), "progress")
dir.create(progress_dir, showWarnings = FALSE)
start_time <- Sys.time()

rep_results <- mclapply(
  seq_len(n_sim),
  run_rep,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

failed <- vapply(rep_results, inherits, logical(1), what = "try-error")
if (any(failed)) {
  warning("worker crashed in reps: ", paste(which(failed), collapse = ", "))
}
rep_results <- rep_results[!failed]

results <- list(
  summary = do.call(rbind, lapply(rep_results, `[[`, "summary")),
  splits = do.call(rbind, lapply(rep_results, `[[`, "splits")),
  nodes = do.call(rbind, lapply(rep_results, `[[`, "nodes")),
  true_splits = true_splits,
  settings = list(
    n_sim = n_sim, n_obs = n_obs, censoring_rate = censoring_rate,
    alpha = alpha, maxdepth = maxdepth, cut = cut, b1 = b1, b2 = b2,
    lambda_c = lambda_c, fit_formula = fit_formula
  )
)

results_dir <- "simulation_custom/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(results, "complex_tree.rds"))
