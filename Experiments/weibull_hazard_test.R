###################################################################
# Pemtree for Weibull distributed survival times ------------------
###################################################################

## setup
source("setup.R")
source("R/sim_weibull.R")
source("R/pemtree_splines.R")

## Helper functions -----------------------------------------------
predict_hazard <- function(object, newdata, ...) {
  newdata$offset <- 0
  as.numeric(
    predict(object$tree, newdata, type = "response")
  )
}

plot_est_hazard <- function(t_grid, hazard_true, models) {
  newdata <- data.frame(tend = t_grid)
  hazard_estimates <- lapply(models, function(model) suppressWarnings(predict_hazard(model, newdata)))

  plot(
    t_grid, hazard_true,
    type = "l",
    col = "black",
    lwd = 3,
    ylim = range(c(hazard_true, unlist(hazard_estimates))),
    xlab = "Time",
    ylab = "Hazard",
    main = "Estimated Hazards from PEM-Trees"
  )
  for (i in seq_along(hazard_estimates)) {
    lines(
      t_grid,
      hazard_estimates[[i]],
      col = adjustcolor("steelblue", alpha.f = 0.3),
      lwd = 1.5
    )
  }
  legend(
    "topleft",
    legend = c("True Weibull hazard", "Estimated hazards"),
    col = c("black", "steelblue"),
    lwd = c(3, 1.5),
    bty = "n"
  )
}

## 3 Scenarios ----------------------------------------------------
# 1) Null scenario: Hazard without level shift -> No split.
# 2) Power scenario: Test for delta = 0.4 and delta = 0.8 if a split is detected.
# 3) Time-varying covariate effect. The effect of a covariate changes at some time point tau.
#    Can we recover a split point here?

## Scenario 1 ------------------------------------------------------
# Parameters for simulation
n <- 1000

shape <- 1.5
rate <- 0.05

censoring_rate <- 0.2
admin_time <- 10

# Repeat model fit 20 times
n_sim <- 20
seeds <- 6531L + seq_len(n_sim)
split_detected <- logical(n_sim)
models <- list()

for (i in seq_len(n_sim)) {
  cat("\n\nRunning simulation: ", i, "\n")
  start_time <- proc.time()[["elapsed"]]
  set.seed(seeds[i])
  sim_data <- sim_weibull(
    n = n,
    shape = shape,
    rate = rate,
    censoring_rate = censoring_rate,
    admin_time = admin_time
  )
  cat("Observed censoring: ", mean(sim_data$status == 0), "\n")

  tree_fit <- pemtree(
    formula = Surv(time, status) ~ bs(tend, df = 5) + offset(offset) | tend,
    data = sim_data,
    min_events = 20,
    maxdepth = 2,
    alpha = 0.05,
    cut = seq(0, admin_time, by = 0.1)
  )

  n_terminal_nodes <- length(nodeids(tree_fit$tree, terminal = TRUE))
  split_detected[i] <- n_terminal_nodes > 1
  models[[length(models) + 1]] <- tree_fit
  run_time <- proc.time()[["elapsed"]] - start_time
  cat("Finished simulation: ", i, ", run time: ", run_time)
}

alpha_hat <- mean(split_detected)
cat("\nEstimated alpha: ", alpha_hat, "\n")

# Plot models
t_grid <- seq(0, 10, length.out = 500)
hazard_true <- rate * shape * t_grid^(shape - 1)
plot_est_hazard(t_grid, hazard_true, models)

## Scenario 2 ------------------------------------------------------
n <- 1000

shape <- 1.5
rate <- 0.05

censoring_rate <- 0.2
admin_time <- 10

tau <- 5
delta <- 0.8

n_sim <- 20
seeds <- 1058L + seq_len(n_sim)

split_detected_power <- logical(n_sim)
models_power <- list()

for (i in seq_len(n_sim)) {
  cat("\n\nRunning simulation: ", i, "\n")
  set.seed(seeds[i])
  start_time <- proc.time()[["elapsed"]]
  sim_data <- sim_weibull(
    n = n,
    shape = shape,
    rate = rate,
    censoring_rate = censoring_rate,
    tau = tau,
    delta = delta,
    admin_time = admin_time
  )
  cat("Observed censoring: ", mean(sim_data$status == 0), "\n")
  tree_fit <- pemtree(
    formula = Surv(time, status) ~ bs(tend, df = 5) + offset(offset) | tend,
    data = sim_data,
    min_events = 20,
    maxdepth = 2,
    alpha = 0.05,
    cut = seq(0, admin_time, by = 0.1)
  )

  n_terminal_nodes <- length(nodeids(tree_fit$tree, terminal = TRUE))
  split_detected_power[i] <- n_terminal_nodes > 1
  models_power[[length(models_power) + 1]] <- tree_fit

  run_time <- proc.time()[["elapsed"]] - start_time
  cat("Finished simulation: ", i, ", run time: ", run_time)
}

power_hat <- mean(split_detected_power)
cat("Estimated Power: ", power_hat, "\n")

# Plot estimated hazards
t_grid <- seq(0, 10, length.out = 500)
hazard_true <- ifelse(
  t_grid <= 5,
  rate * shape * t_grid^(shape - 1),
  rate * shape * t_grid^(shape - 1) * exp(delta)
)
plot_est_hazard(t_grid, hazard_true, models_power)

## Scenario 3 ------------------------------------------------------
