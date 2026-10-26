### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/permtree.R")

## Helper functions ---------------------------------------------------
run_simulation <- function(n_sim, n_obs, n_pilot,
                           seeds, calibration_seed,
                           formula, covariate_spec, tau = NULL,
                           admin_time, censoring_rate) {
  set.seed(calibration_seed)
  breaks <- if (is.null(tau)) numeric(0) else tau
  cat("\n\nRunning calibration")
  covar <- lapply(covariate_spec, function(s) {
    assert_true("distfun" %in% names(s))
    fun <- s$distfun
    args <- s[names(s) != "distfun"]
    do.call(fun, c(list(n = n_pilot), args))
  })
  data_pilot <- as.data.frame(covar)
  data_pilot$id <- seq_len(n_pilot)
  event_times <- sim_event_times(formula, data_pilot, breaks = breaks)
  lambda_c <- calibrate_censoring(event_times, censoring_rate, admin_time)

  fit_cut <- sort(unique(c(
    seq(0, admin_time, by = 0.1), tau
  )))
  fit_formula <- Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | z

  split_detected_perm <- logical(n_sim)
  split_detected_pem <- logical(n_sim)

  p_vals_perm <- numeric(n_sim)
  p_vals_pem <- numeric(n_sim)

  coef_perm <- list()
  coef_pem <- list()

  for (i in seq_len(n_sim)) {
    cat("\n\nRunning simulation: ", i, "with seed: ", seeds[i], "\n")
    set.seed(seeds[i])
    start_time <- proc.time()[["elapsed"]]

    data <- as.data.frame(lapply(covariate_spec, function(s) {
      fun <- s$distfun
      args <- s[names(s) != "distfun"]
      do.call(fun, c(list(n = n_obs), args))
    }))
    data$id <- seq_len(n_obs)
    sim_data <- sim_surv(
      formula = formula,
      data = data,
      breaks = breaks,
      admin_time = admin_time,
      lambda_c = lambda_c
    )
    sim_data$z <- factor(sim_data$z, levels = c(0, 1))
    cat("Observed censoring rate: ", mean(sim_data$status == 0), "\n")

    # Fit trees
    pemtree_fit <- pemtree(
      formula = fit_formula,
      data = sim_data,
      cut = fit_cut,
      maxdepth = 2,
      alpha = 0.05
    )

    permtree_fit <- PErMtree(
      formula = fit_formula,
      data = sim_data,
      n_perm = 199L,
      alpha = 0.05,
      cut = fit_cut
    )

    split_detected_pem[i] <- length(nodeids(pemtree_fit$tree, terminal = TRUE)) > 1
    split_detected_perm[i] <- permtree_fit$terminal_nodes > 1

    root <- node_party(pemtree_fit$tree)
    p_vals_pem[i] <- info_node(root)$test["p.value", ]
    p_vals_perm[i] <- permtree_fit$p_vals[["z"]]

    coef_pem[[i]] <- coef(pemtree_fit$tree)
    coef_perm[[i]] <- lapply(permtree_fit$models, coef)

    runtime <- proc.time()[["elapsed"]] - start_time
    cat("Finished simulation: ", i, "with runtime: ", round(runtime, 4))
  }

  list(
    split_detected = list(permtree = split_detected_perm, pemtree = split_detected_pem),
    alpha_hat = list(permtree = mean(split_detected_perm), pemtree = mean(split_detected_pem)),
    p_val = list(permtree = p_vals_perm, pemtree = p_vals_pem),
    coefficients = list(permtree = coef_perm, pemtree = coef_pem)
  )
}

## Comparison of PErMtree and pemtree ---------------------------------
## 3 scenarios are checked here
## Scenario 1: Null scenario with no interactions
## Scenario 2: Power Scenario with an effect of 0.5
## Scenario 3: More complex scenario with two possible split variables
##             -> One has an effect and one has no effect

## Global parameters ---------------------------------------------------
n_pilot <- 10000
n_obs <- 500
shape <- 1.5
rate <- 0.05

log_hazard_base <- function(t) {
  log(shape) + log(rate) + (shape - 1) * log(t)
}

beta1 <- 0.7
formula <- ~ log_hazard_base(t) + beta1 * x + delta * x * (z == 1)

admin_time <- 10
censoring_rate <- 0.3

n_sim <- 50

## Scenario 1 ---------------------------------------------------------
covariate_spec <- list(
  x = list(distfun = rnorm),
  z = list(distfun = rbinom, size = 1, prob = 0.5)
)
delta <- 0

seeds <- 1136 + seq_len(n_sim)
calibration_seed <- 7

results_1 <- run_simulation(
  n_sim = n_sim, n_obs = n, n_pilot = n_pilot,
  seeds = seeds, calibration_seed = calibration_seed,
  formula = formula, covariate_spec = covariate_spec,
  tau = NULL, admin_time = admin_time, censoring_rate = censoring_rate
)

results_1$alpha_hat


## Scenario 2 ---------------------------------------------------------
covariate_spec <- list(
  x = list(distfun = rnorm),
  z = list(distfun = rbinom, size = 1, prob = 0.5)
)
delta <- 0.5

seeds <- 1257 + seq_len(n_sim)
calibration_seed <- 8

results_2 <- run_simulation(
  n_sim = n_sim, n_obs = n, n_pilot = n_pilot,
  seeds = seeds, calibration_seed = calibration_seed,
  formula = formula, covariate_spec = covariate_spec,
  tau = NULL, admin_time = admin_time, censoring_rate = censoring_rate
)

results_2$alpha_hat


## Scenario 3 ---------------------------------------------------------
covariate_spec <- list(
  x = list(distfun = rnorm),
  z1 = list(distfun = rbinom, size = 1, prob = 0.6),
  z2 = list(distfun = rbinom, size = 1, prob = 0.3)
)
delta1 <- 0.4
delta2 <- 0.0

formula <- ~ log_hazard_base(t) + beta1 * x + delta1 * x * (z1 == 1) + delta2 * x * (z2 == 1)

seeds <- 20 + seq_len(n_sim)
calibration_seed <- 9

split_detected <- list(permtree = logical(n_sim), pemtree = logical(n_sim))
split_vars <- list(permtree = character(n_sim), pemtree = character(n_sim))
p_vals_adjusted <- list(permtree = list(), pemtree = list())


# Calibrate once for this scenario
set.seed(calibration_seed)
cat("\n\nRunning calibration")
covar <- lapply(covariate_spec, function(s) {
  assert_true("distfun" %in% names(s))
  fun <- s$distfun
  args <- s[names(s) != "distfun"]
  do.call(fun, c(list(n = n_pilot), args))
})
data_pilot <- as.data.frame(covar)
data_pilot$id <- seq_len(n_pilot)
event_times <- sim_event_times(formula, data_pilot)
lambda_c <- calibrate_censoring(event_times, censoring_rate, admin_time)

# Model specification
fit_cut <- seq(0, admin_time, by = 0.1)
fit_formula <- Surv(time, status) ~ bs(tend, df = 6, degree = 3) + x + offset(offset) |
  z1 + z2

for (i in seq_len(n_sim)) {
  cat("\nRunning simulation", i, "with seed", seeds[i], "\n")
  set.seed(seeds[i])
  start_time <- proc.time()[["elapsed"]]

  covar <- lapply(covariate_spec, function(s) {
    fun <- s$distfun
    args <- s[names(s) != "distfun"]
    do.call(fun, c(list(n = n), args))
  })
  sim_input <- as.data.frame(covar)
  sim_input$id <- seq_len(n)

  sim_data <- sim_surv(
    formula = formula,
    data = sim_input,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
  sim_data$z1 <- factor(sim_data$z1, levels = c(0, 1))
  sim_data$z2 <- factor(sim_data$z2, levels = c(0, 1))
  cat("Observed censoring rate: ", mean(sim_data$status == 0), "\n")

  # Fit both methods
  permtree_fit <- PErMtree(
    formula = fit_formula,
    data = sim_data,
    n_perm = 199L,
    alpha = 0.05,
    cut = fit_cut
  )
  pemtree_fit <- pemtree(
    formula = fit_formula,
    data = sim_data,
    cut = fit_cut,
    maxdepth = 2,
    alpha = 0.05
  )

  perm_split <- permtree_fit$split
  tree_pem <- pemtree_fit$tree
  root_pem <- node_party(tree_pem)
  pem_split <- length(nodeids(tree_pem, terminal = TRUE)) > 1
  split_detected$permtree[i] <- perm_split
  split_detected$pemtree[i] <- pem_split

  split_vars$permtree[i] <- if (perm_split) permtree_fit$split_var else NA_character_
  if (pem_split) {
    var_id <- root_pem$split$varid
    split_vars$pemtree[i] <- colnames(data_party(tree_pem))[var_id]
  } else {
    split_vars$pemtree[i] <- NA_character_
  }

  p_vals_adjusted$permtree[[i]] <- permtree_fit$p_adjusted
  p_vals_pem <- info_node(root_pem)$test["p.value", , drop = TRUE]
  p_vals_adjusted$pemtree[[i]] <- p.adjust(p_vals_pem, method = "bonferroni")

  run_time <- proc.time()[["elapsed"]] - start_time
  cat("Finished simulation: ", i, "with runtime: ", round(run_time, 2))
}
