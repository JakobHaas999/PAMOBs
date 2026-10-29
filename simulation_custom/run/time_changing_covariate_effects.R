### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
results_dir <- "simulation_custom/results"

## Helper functions ---------------------------------
strip_tree <- function(tree) {
  nodes <- lapply(as.list(node_party(tree)), function(n) {
    n$info$object <- NULL
    n$info$estfun <- NULL
    n
  })
  tree$node <- as.partynode(nodes)
  tree$info$call <- NULL
  tree$info$fit <- NULL
  tree
}

predict_hazard <- function(tree, newdata) {
  newdata$offset <- 0
  newdata$node <- predict(tree, newdata = newdata, type = "node")
  models <- nodeapply(tree,
    ids = unique(newdata$node),
    FUN = function(n) n$info$object
  )
  newdata$log_hazard <- NA_real_
  for (id in names(models)) {
    rows <- newdata$node == as.integer(id)
    newdata$log_hazard[rows] <- suppressWarnings(
      predict(models[[id]], newdata = newdata[rows, , drop = FALSE], type = "link")
    )
  }
  newdata$hazard <- exp(newdata$log_hazard)
  newdata
}

## Global parameters --------------------------------

n_sim <- 100
n_obs <- 1000
n_pilot <- 10000

shape <- 1.5
rate <- 0.05

censoring_rate <- 0.3
admin_time <- 10
fit_cut <- seq(0, admin_time, by = 0.1)

log_hazard0 <- function(t) log(shape) + log(rate) + (shape - 1) * log(t)

## Scenario 1 ----------------------------------------
## Simple covariate changing effect
beta1 <- 0.4
beta2 <- 0.6
delta <- 0.8
tau <- 4

formula_1 <- ~ log_hazard0(t) + beta1 * x1 + beta2 * x2 + delta * x1 * (t >= tau)

# Calibrate censoring lambda
set.seed(1)
data_pilot_1 <- data.frame(
  x1 = rnorm(n_pilot),
  x2 = runif(n_pilot, -1, 1),
  id = seq_len(n_pilot)
)
event_times_1 <- sim_event_times(
  formula = formula_1,
  data = data_pilot_1,
  breaks = tau
)
lambda_c_1 <- calibrate_censoring(event_times_1, censoring_rate, admin_time)

seeds_1 <- 1309 + seq_len(n_sim)
nd <- expand.grid(tend = fit_cut[-1], x1 = c(-1, 0, 1), x2 = 0)

run_one_1 <- function(i) {
  message("Simulation ", i, " of ", n_sim, " started")
  set.seed(seeds_1[i])

  data <- data_frame(
    id = seq_len(n_obs),
    x1 = rnorm(n_obs),
    x2 = runif(n_obs, -1, 1)
  )
  sim_data <- sim_surv(
    formula = formula_1,
    data = data,
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c_1
  )

  tree_fit <- pemtree(
    formula = Surv(time, status) ~ bs(tend, df = 6) + x1 + x2 + offset(offset) | tend,
    data = sim_data,
    cut = fit_cut,
    maxdepth = 2
  )

  list(
    hazards = cbind.data.frame(sim = i, predict_hazard(tree_fit$tree, nd)),
    tree = strip_tree(tree_fit$tree),
    censoring = mean(sim_data$status == 0)
  )
}

out_1 <- parallel::mclapply(seq_len(n_sim), run_one_1, mc.cores = 2)
failed <- vapply(out_1, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Simulations failed: ", paste(which(failed), collapse = ", "))



dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(trees = trees_1, hazards = hazards_1, seeds = seeds_1),
  file.path(results_dir, "time_changing_effects_1.rds")
)

## Scenario 1 ----------------------------------------
## Two covariate changing effects
