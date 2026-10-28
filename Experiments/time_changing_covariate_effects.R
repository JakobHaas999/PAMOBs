##########################################################
## Time-changing covariate effects for pamtrees ##########
##########################################################

## setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

## 3 scenarios --------------------------------------------
## Scenario 1: simple covariate changing effect
## Scenario 2: Two covariate changing effects
## Scenario 3: Two covariates have an changing effect


## Global parameters
n_sim <- 50
n_obs <- 1000
n_pilot <- 5000

shape <- 1.5
rate <- 0.05

censoring_rate <- 0.3
admin_time <- 10


## Scenario 1 --------------------------------------------

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -1, max = 1)
)

log_hazard0 <- function(t) {
  log(shape) + log(rate) + (shape - 1) * log(t)
}

beta1 <- 0.4
beta2 <- 0.6
delta <- 0.8
tau <- 4

formula_1 <- ~ log_hazard0(t) + beta1 * x1 + beta2 * x2 + delta * x1 * (t >= tau)

# Calibrate censoring lambda
set.seed(1)
data_pilot <- as.data.frame(lapply(covariate_spec, function(s) {
  distfun <- s$distfun
  args <- s[names(s) != "distfun"]
  do.call(distfun, c(list(n = n_pilot), args))
}))
data_pilot$id <- seq_len(n_pilot)
event_times <- sim_event_times(
  formula = formula_1,
  data = data_pilot,
  breaks = tau
)
lambda_c <- calibrate_censoring(event_times, censoring_rate, admin_time)

split_detected <- logical(n_sim)
tau_hats <- numeric(n_sim)
trees <- list()

seeds <- 2147 + seq_len(n_sim)

for (i in seq_len(n_sim)) {
  cat("\n\nRunning simulation : ", i, "with seed: ", seeds[i], "\n")
  set.seed(seeds[i])

  data <- as.data.frame(lapply(covariate_spec, function(s) {
    distfun <- s$distfun
    args <- s[names(s) != "distfun"]
    do.call(distfun, c(list(n = n_obs), args))
  }))
  data$id <- seq_len(n_obs)
  sim_data <- sim_surv(
    formula = formula_1,
    data = data,
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
  cat("Observed censoring rate: ", mean(sim_data$status == 0), "\n")

  fit_cut <- seq(0, admin_time, by = 0.1)
  tree_fit <- pemtree(
    formula = Surv(time, status) ~ bs(tend, df = 6) + x1 + x2 + offset(offset) | tend,
    data = sim_data,
    cut = fit_cut,
    maxdepth = 2
  )

  tree <- tree_fit$tree
  n_terminal_nodes <- length(nodeids(tree, terminal = TRUE))
  split_detected[[i]] <- n_terminal_nodes > 1
  tau_hats[[i]] <- if (n_terminal_nodes > 1) {
    breaks_split(split_node(node_party(tree)))
  } else {
    NA_real_
  }
  trees[[i]] <- tree
}

mean(split_detected)
summary(tau_hats)

## Scenario 2 --------------------------------------------

log_hazard0 <- function(t) {
  log(shape) + log(rate) + (shape - 1) * log(t)
}

tau1 <- 3
tau2 <- 5
breaks_2 <- c(tau1, tau2)

delta1 <- 0.6
delta2 <- 0.6

beta1 <- 0.4
beta2 <- 0.6

formula_2 <- ~ log_hazard0(t) + beta1 * x1 + beta2 * x2 + delta1 * x1 * (t >= tau1) + delta2 * x2 * (t >= tau2)

# Calibration
set.seed(2)
data_pilot <- data.frame(
  id = seq_len(n_pilot),
  x1 = rnorm(n_pilot),
  x2 = runif(n_pilot, -1, 1)
)
event_times <- sim_event_times(formula_2, data_pilot, breaks_2)
lambda_c <- calibrate_censoring(event_times, censoring_rate, admin_time)

seeds_2 <- 24 + seq_len(n_sim)

tau_hats_2 <- vector("list", n_sim)
n_splits_2 <- integer(n_sim)
trees_2 <- vector("list", n_sim)

fit_cut <- seq(0, admin_time, by = 0.1)

for (i in seq_len(n_sim)) {
  cat("\nRunning simulation:", i, "with seed", seeds_2[i], "\n")
  set.seed(seeds_2[i])

  data <- data.frame(
    id = seq_len(n_obs),
    x1 = rnorm(n_obs),
    x2 = runif(n_obs, -1, 1)
  )
  sim_data <- sim_surv(
    formula = formula_2,
    data = data,
    breaks = breaks_2,
    admin_time = admin_time,
    lambda_c = lambda_c
  )

  cat("Observed censoring:", mean(sim_data$status == 0), "\n")

  tree_fit <- pemtree(
    formula = Surv(time, status) ~ bs(tend, df = 6) + x1 + x2 + offset(offset) | tend,
    data = sim_data,
    cut = fit_cut,
    maxdepth = 3,
    alpha = 0.05
  )

  tree <- tree_fit$tree

  inner_ids <- setdiff(nodeids(tree), nodeids(tree, terminal = TRUE))
  tau_hats_2[[i]] <- if (length(inner_ids)) {
    sort(unlist(nodeapply(
      tree,
      ids = inner_ids,
      FUN = function(node) {
        breaks_split(split_node(node))
      }
    ), use.names = FALSE))
  } else {
    numeric(0)
  }

  n_splits_2[[i]] <- length(tau_hats_2[[i]])
  trees_2[[i]] <- tree
}

table(n_splits_2)
