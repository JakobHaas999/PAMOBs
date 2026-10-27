############################################################
## Analysis of pamtree #####################################

### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pamtree.R")

## Simulation 1: No level shift

# Simulate data
n_sim <- 20

rate <- 0.05
shape <- 1.5
log_hazard <- function(t) {
  log(rate) + log(shape) + (shape - 1) * log(t)
}
formula <- ~ log_hazard(t)

n_pilot <- 50000
n_obs <- 1000

admin_time <- 10
censoring_rate <- 0.3


# Calibrate censoring lambda
set.seed(101)
data_pilot <- data.frame(id = seq_len(n_pilot))
event_times <- sim_event_times(formula, data_pilot)
lambda_c <- calibrate_censoring(event_times, censoring_rate, admin_time)
lambda_c

# Calculate optimal smoothing parameter
set.seed(102)
data_sp <- sim_surv(formula, data.frame(id = seq_len(n_obs)), admin_time = admin_time, lambda_c = lambda_c)

cut <- seq(0, admin_time, by = 0.1)

set.seed(103)
grid_search_res <- pam_grid_search(
  formula = Surv(time, status) ~ s(tend, bs = "ps") + offset(offset),
  data = data_sp,
  cut = cut,
  sp_space = c(0, 0.1, 1, 10, 100, 1000),
  method = "cv",
  folds = 3
)

sp_fixed <- grid_search_res$selected_sp
sp_fixed

# Run simulation
seeds <- 1511L + seq_len(n_sim)
split_detected <- logical(n_sim)
trees <- list()
wrong_taus <- c()

for (i in seq_len(n_sim)) {
  cat("\nRunning simulation: ", i, "with seed: ", seeds[i], "\n")
  set.seed(seeds[i])

  sim_data <- sim_surv(
    formula = formula,
    data = data.frame(id = seq_len(n_obs)),
    admin_time = admin_time,
    lambda_c = lambda_c
  )
  cat("Observed censoring rate: ", mean(sim_data$status == 0), "\n")

  tree_fit <- pamtree(
    formula = Surv(time, status) ~ s(tend, bs = "ps") + offset(offset) | tend,
    data = sim_data,
    cut = cut,
    sp = sp_fixed,
    maxdepth = 2,
    alpha = 0.05
  )

  n_terminal_nodes <- length(nodeids(tree_fit$tree, terminal = TRUE))
  split_detected[i] <- n_terminal_nodes > 1
  trees[[i]] <- tree_fit$tree
  if (n_terminal_nodes > 1) {
    root <- node_party(tree_fit$tree)
    wrong_taus <- c(wrong_taus, breaks_split(split_node(root)))
  }
}

mean(split_detected)
mean(wrong_taus)

## Simulation 2: no level shift with a constant baseline hazard
