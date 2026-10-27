##############################################################################
## Level shift scenarios for hazard rates of different distributions #########
##############################################################################

## setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

## Helper functions ----------------------------------------------------------
run_simulation <- function(n_sim, n_obs, n_pilot,
                           seeds, calibration_seed, formula,
                           tau, admin_time, censoring_rate, interval, df = 6) {
  set.seed(calibration_seed)
  event_times <- sim_event_times(formula,
    data = data.frame(id = seq_len(n_pilot)),
    breaks = tau, interval = interval
  )
  lambda_c <- calibrate_censoring(event_times, censoring_rate, admin_time)

  split_detected <- logical(n_sim)
  trees <- vector("list", n_sim)
  tau_hats <- numeric(n_sim)

  for (i in seq_len(n_sim)) {
    cat("\n\nRunning simulation: ", i, "with seed: ", seeds[i], "\n")
    set.seed(seeds[i])
    start_time <- proc.time()[["elapsed"]]

    sim_data <- sim_surv(
      formula = formula,
      data = data.frame(id = seq_len(n_obs)),
      breaks = tau,
      admin_time = admin_time,
      lambda_c = lambda_c,
      interval = interval
    )
    cat("Observed censoring rate: ", mean(sim_data$status == 0), "\n")

    fit_cut <- sort(unique(c(seq(0, admin_time, by = 0.1), tau)))
    tree_fit <- pemtree(
      formula = Surv(time, status) ~ bs(tend, df = df) + offset(offset) | tend,
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
      tau_hat <- breaks_split(split_node(root))
    } else {
      tau_hat <- NA_real_
    }
    tau_hats[i] <- tau_hat

    run_time <- proc.time()[["elapsed"]] - start_time
    cat("Finished simulation: ", i, "with runtime: ", round(run_time, 4))
  }

  list(
    split_detected = split_detected,
    tau_hats = tau_hats,
    trees = trees,
    power = mean(split_detected)
  )
}


## 3 scenarios with level shifts ---------------------------------------------
## Scenario 1: Cubic baseline hazard
## Scenario 2: Weibull baseline hazard
## Scenario 3: Log-Normal baseline hazard

## Scenario 1 ----------------------------------------------------------------
coeff <- c(-2.8, 0.25, -0.025, 0.001)
delta <- log(3)
tau <- 5
log_hazard <- function(t) {
  coeff[1] + coeff[2] * t + coeff[3] * t^2 + coeff[4] * t^3 + delta * (t >= tau)
}

formula <- ~ log_hazard(t)

n_sim <- 50
n_obs <- 1000
n_pilot <- 10000...............

seeds <- 8L + seq_len(n_sim)
calibration_seed <- 4

admin_time <- 10
censoring_rate <- 0.3
interval <- c(1e-8, 30)

results_1 <- run_simulation(
  n_sim, n_obs, n_pilot,
  seeds, calibration_seed,
  formula, tau, admin_time,
  censoring_rate, interval
)

power_1 <- results_1$power
cat("Scenario 1 power: ", power_1, "\n")

summary(results_1$tau_hats)
hist(
  results_1$tau_hats,
  breaks = 20,
  main = "Estimated change points - Cubic hazard",
  xlab = expression(hat(tau))
)
abline(v = tau, col = "red", lwd = 2)

## Bias and RMSE
tau_hat <- na.omit(results_1$tau_hats)
bias_tau_1 <- mean(tau_hat - tau)
rmse_tau_1 <- sqrt(mean((tau_hat - tau)^2))
cat("Scenario 1 Bias:", bias_tau_1, "\n")
cat("Scenario 1 RMSE:", rmse_tau_1, "\n")


## Scenario 2 ----------------------------------------------------------------
shape <- 1.5
rate <- 0.05
delta <- log(3)
tau <- 5

log_hazard <- function(t) {
  log(rate) + log(shape) + (shape - 1) * log(t) + delta * (t >= tau)
}
formula <- ~ log_hazard(t)

n_sim <- 50
n_obs <- 1000
n_pilot <- 10000

seeds <- 1626L + seq_len(n_sim)
calibration_seed <- 5

admin_time <- 10
censoring_rate <- 0.3
interval <- c(1e-8, 30)

results_2 <- run_simulation(
  n_sim, n_obs, n_pilot,
  seeds, calibration_seed,
  formula, tau, admin_time,
  censoring_rate, interval
)

power_2 <- results_2$power
cat("Scenario 2 power: ", power_2, "\n")

summary(results_2$tau_hats)
hist(
  results_2$tau_hats,
  breaks = 20,
  main = "Estimated change points - Weibull hazard",
  xlab = expression(hat(tau))
)
abline(v = tau, col = "red", lwd = 2)

## Bias and RMSE
tau_hat <- na.omit(results_2$tau_hats)
bias_tau_2 <- mean(tau_hat - tau)
rmse_tau_2 <- sqrt(mean((tau_hat - tau)^2))
cat("Scenario 2 Bias:", bias_tau_2, "\n")
cat("Scenario 2 RMSE:", rmse_tau_2, "\n")


## Scenario 3 ----------------------------------------------------------------
meanlog <- log(6)
sdlog <- 0.7
tau <- 5
delta <- log(3)

log_hazard <- function(t) {
  dlnorm(t, meanlog = meanlog, sdlog = sdlog, log = TRUE) -
    plnorm(t, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE, log.p = TRUE) +
    delta * (t >= tau)
}

formula <- ~ log_hazard(t)

n_sim <- 50
n <- 1000
n_pilot <- 10000

seeds <- 1903L + seq_len(n_sim)
calibration_seed <- 6

admin_time <- 10
censoring_rate <- 0.3
interval <- c(1e-8, 30)

results_3 <- run_simulation(
  n_sim, n, n_pilot, seeds, calibration_seed,
  formula, tau, admin_time, censoring_rate, interval,
  df = 5
)

power_3 <- results_3$power
cat("Scenario 3 power: ", power_3, "\n")

summary(results_3$tau_hats)
hist(
  results_3$tau_hats,
  breaks = 20,
  main = "Estimated change points - Lognormal hazard",
  xlab = expression(hat(tau))
)
abline(v = tau, col = "red", lwd = 2)

## Bias and RMSE
tau_hat <- na.omit(results_3$tau_hats)
bias_tau_3 <- mean(tau_hat - tau)
rmse_tau_3 <- sqrt(mean((tau_hat - tau)^2))
cat("Scenario 3 Bias:", bias_tau_3, "\n")
cat("Scenario 3 RMSE:", rmse_tau_3, "\n")



