##############################################################################
## Null scenarios for hazard rates of different distributions ################
##############################################################################

## setup ---------------------------------------------------------------------

source("setup.R")
source("R/pemtree_splines.R")
source("R/sim_surv.R")

## Helper functions ----------------------------------------------------------
predict_hazard <- function(fit, newdata) {
  fit <- if (inherits(fit, "pemtree")) {
    fit$tree
  } else {
    fit
  }
  newdata$offset <- 0
  as.numeric(
    predict(fit, newdata = newdata, type = "response")
  )
}

plot_hazard_est <- function(t_grid, models, true_hazard) {
  newdata <- data.frame(tend = t_grid)
  hazard_estimates <- lapply(models, function(model) {
    suppressWarnings(predict_hazard(model, newdata))
  })

  plot(
    t_grid, true_hazard,
    type = "l",
    col = "black",
    lwd = 3,
    ylim = range(c(true_hazard, unlist(hazard_estimates))),
    xlab = "t", ylab = expression(lambda(t))
  )
  for (i in seq_along(hazard_estimates)) {
    lines(
      t_grid, hazard_estimates[[i]],
      col = adjustcolor("steelblue", alpha.f = 0.3),
      lwd = 1.5
    )
  }
}

calc_ise <- function(t_grid, hazard_est, hazard_true) {
  errors_sq <- (hazard_est - hazard_true)^2

  sum(diff(t_grid) * (head(errors_sq, -1) + tail(errors_sq, -1)) / 2)
}

run_simulation <- function(n_sim, n_obs, n_pilot,
                           seeds, calibration_seed,
                           formula, admin_time, censoring_rate) {
  set.seed(calibration_seed)
  event_times <- sim_event_times(formula, data.frame(id = seq_len(n_pilot)))
  lambda_c <- calibrate_censoring(event_times, censoring_rate = censoring_rate, admin_time = admin_time)

  split_detected <- logical(n_sim)
  trees <- list()
  wrong_taus <- c()

  for (i in seq_len(n_sim)) {
    cat("\n\nRunning simulation: ", i, "with seed: ", seeds[i], "\n")
    set.seed(seeds[i])
    start_time <- proc.time()[["elapsed"]]

    sim_data <- sim_surv(
      formula = formula,
      data = data.frame(id = seq_len(n)),
      admin_time = admin_time,
      lambda_c = lambda_c
    )
    cat("Observed censoring rate: ", mean(sim_data$status == 0), "\n")

    fit_cut <- seq(0, admin_time, by = 0.1)
    tree_fit <- pemtree(
      formula = Surv(time, status) ~ bs(tend, df = 6) + offset(offset) | tend,
      data = sim_data,
      cut = fit_cut,
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
    runtime <- proc.time()[["elapsed"]] - start_time
    cat("Finished simulation: ", i, "with runtime: ", round(runtime, 4))
  }

  list(
    split_detected = split_detected, alpha_hat = mean(split_detected),
    trees = trees, wrong_taus = wrong_taus
  )
}

## 3 scenarios ---------------------------------------------------------------
## Scenario 1: Cubic baseline hazard
## Scenario 2: Weibull baseline hazard
## Scenario 3: Log-Normal baseline hazard

## Scenario 1 ----------------------------------------------------------------
coeff <- c(-2.8, 0.25, -0.025, 0.001)
log_hazard <- function(t) coeff[1] + coeff[2] * t + coeff[3] * t^2 + coeff[4] * t^3
formula <- ~ log_hazard(t)
n <- 1000
admin_time <- 10
censoring_rate <- 0.3
n_pilot <- 10000
n_sim <- 50
seeds <- 1536 + seq_len(n_sim)

results_1 <- run_simulation(
  n_sim = n_sim, n_obs = n, n_pilot = n_pilot,
  seeds = seeds, calibration_seed = 1,
  formula = formula, admin_time = admin_time,
  censoring_rate = censoring_rate
)

alpha_hat_1 <- results_1$alpha_hat
cat("Scenario 1 alpha hat: ", alpha_hat_1, "\n")

# Plot predictions
t_grid <- seq(0.01, ceiling(max(sim_data$time)), length.out = 500)
true_hazard_1 <- exp(log_hazard(t_grid))
plot_est_hazard(t_grid, true_hazard_1, results_1$trees)

# Calculate ISE
hazard_est_1 <- lapply(results_1$trees, function(m) {
  suppressWarnings(predict_hazard(m, newdata = data.frame(tend = t_grid)))
})
log_ise_1 <- vapply(hazard_est_1, function(h) {
  calc_ise(t_grid, log(h), log(true_hazard_1))
}, FUN.VALUE = numeric(1))
mean(log_ise_1)

results_1$wrong_taus


## Scenario 2 ----------------------------------------------------------------
shape <- 1.5
rate <- 0.05
formula <- ~ log(rate) + log(shape) + (shape - 1) * log(t)
n <- 1000
admin_time <- 10
censoring_rate <- 0.3
n_pilot <- 10000
seeds <- 1850L + seq_len(n_sim)

results_2 <- run_simulation(
  n_sim = n_sim, n_obs = n, n_pilot = n_pilot,
  seeds = seeds, calibration_seed = 2, formula = formula,
  admin_time = admin_time, censoring_rate = censoring_rate
)

alpha_hat_2 <- results_2$alpha_hat
cat("Scenario 2 alpha hat: ", alpha_hat_2, "\n")

# Plot predictions
t_grid <- seq(0.01, ceiling(max(sim_data$time)), length.out = 500)
true_hazard_2 <- rate * shape * t_grid^(shape - 1)
plot_est_hazard(t_grid, true_hazard_2, results_2$trees)

# Calculate ISE
hazard_est_2 <- lapply(results_2$trees, function(m) {
  suppressWarnings(predict_hazard(m, newdata = data.frame(tend = t_grid)))
})
log_ise_2 <- vapply(hazard_est_2, function(h) {
  calc_ise(t_grid, log(h), log(true_hazard_2))
}, FUN.VALUE = numeric(1))
mean(log_ise_2)

results_2$wrong_taus

## Scenario 3
