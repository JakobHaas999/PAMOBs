### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

### 2 Scenarios
### Scenario 1: Different weibull parameter shape between z = 0 and z = 1
### Scenario 2: Different weibull parameter shape between numerical z <= c and z > c

### Helper functions


### Global parameters
shape0 <- 1.5
rate0 <- 0.05
shape1 <- 2.5
rate1 <- 0.05^(shape1 / 1.5)

# Plot hazard
tvals <- seq(0, 10, length.out = 500)
weibull_hazard <- function(t, shape, rate) exp(log(rate) + log(shape) + (shape - 1) * log(t))
plot(tvals, weibull_hazard(tvals, shape0, rate0),
  type = "l",
  ylim = c(0, weibull_hazard(10, shape1, rate1)),
  xlab = "t", ylab = expression(lambda(t)), lwd = 2
)
lines(tvals, weibull_hazard(tvals, shape1, rate1), col = "red", lwd = 2)
grid()

beta <- 0.6

admin_time <- 10
censoring_rate <- 0.4

n_sim <- 20
n_pilot <- 5000
n_obs <- 1000

fit_cut <- seq(0, admin_time, by = 0.1)
fit_formula <- Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | z

# Scenario 1 -------------------------------------------------------------
log_hazard <- function(t, z) {
  h0 <- log(shape0) + log(rate0) + (shape0 - 1) * log(t)
  h1 <- log(shape1) + log(rate1) + (shape1 - 1) * log(t)
  ifelse(z == 1, h1, h0)
}

sim_formula <- ~ log_hazard(t, z) + beta * x

# Calibration
set.seed(12)
data_pilot <- data.frame(
  id = seq_len(n_pilot),
  x = rnorm(n_pilot),
  z = factor(rbinom(n_pilot, size = 1, prob = 0.5), levels = c(0, 1))
)
event_times_pilot <- sim_event_times(
  formula = sim_formula,
  data = data_pilot,
  maxt = admin_time
)
lambda_c <- calibrate_censoring(event_times_pilot, censoring_rate, admin_time)
cat(sprintf("Censoring distribution lambda in scenario 1: %.4f\n", lambda_c))

trees_1 <- vector("list", length = n_sim)
split_detected_1 <- logical(n_sim)
root_p_value_1 <- numeric(n_sim)
observed_censoring_1 <- numeric(n_sim)
errors_1 <- character()

seeds <- 816 + seq_len(n_sim)

for (i in seq_len(n_sim)) {
  message("Running simulation ", i, " with seed ", seeds[i])
  set.seed(seeds[i])

  data <- data.frame(
    id = seq_len(n_obs),
    x = rnorm(n_obs),
    z = factor(rbinom(n_obs, size = 1, prob = 0.5), levels = c(0, 1))
  )
  sim_data <- sim_surv(
    formula = sim_formula,
    data = data,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
  message("Observed censoring rate: ", mean(sim_data$status == 0))
  observed_censoring_1[i] <- mean(sim_data$status == 0)

  fit <- tryCatch(
    pemtree(
      formula = fit_formula,
      data = sim_data,
      cut = fit_cut,
      maxdepth = 2,
      alpha = 0.05
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    errors_1 <- c(errors_1, conditionMessage(fit))
    trees_1[[i]] <- NA
    split_detected_1[i] <- NA
    root_p_value_1[i] <- NA
    next
  }

  tree <- fit$tree
  trees_1[[i]] <- tree
  split_detected_1[i] <- length(nodeids(tree, terminal = TRUE)) > 1
  root <- node_party(tree)
  root_p_value_1[i] <- info_node(root)$p.value
}

power_1 <- mean(split_detected_1, na.rm = TRUE)
cat(sprintf("Power in scenario 1: %.2f\n", power_1))

# Plot estimated hazards
t_grid <- fit_cut[-1]
newdata <- data.frame(tend = t_grid, x = 0, offset = 0)
hazard_est_1 <- do.call(rbind, lapply(seq_along(trees_1), function(i) {
  tree <- trees_1[[i]]
  if (!inherits(tree, "modelparty")) {
    return(NULL)
  }
  nodes <- predict(tree, newdata = data.frame(z = factor(0:1, levels = c(0, 1))), type = "node")
  do.call(rbind, lapply(1:2, function(j) {
    m <- refit.modelparty(tree, node = nodes[j])
    data.frame(sim = i, z = j - 1, t = t_grid, hazard = predict(m, newdata = newdata, type = "response"))
  }))
}))

p1 <- ggplot(hazard_est_1, aes(x = t, y = hazard, group = interaction(sim, z), colour = factor(z))) +
  geom_step(alpha = 0.3) +
  stat_function(fun = function(t) weibull_hazard(t, shape0, rate0), colour = "black", inherit.aes = FALSE) +
  stat_function(fun = function(t) weibull_hazard(t, shape1, rate1), colour = "red", inherit.aes = FALSE) +
  scale_colour_manual(values = c("0" = "grey40", "1" = "salmon"), name = "z") +
  labs(y = "Hazard (x = 0)") +
  theme_bw()
print(p1)

root_p_value_1

# Scenario 2 -------------------------------------------------------------
log_hazard <- function(t, z) {
  h0 <- log(shape0) + log(rate0) + (shape0 - 1) * log(t)
  h1 <- log(shape1) + log(rate1) + (shape1 - 1) * log(t)
  ifelse(z <= 0, h0, h1)
}

sim_formula <- ~ log_hazard(t, z) + beta * x

# Calibration
set.seed(13)
data_pilot <- data.frame(
  id = seq_len(n_pilot),
  x = rnorm(n_pilot),
  z = round(runif(n_pilot, min = -1, max = 1), 2)
)
event_times_pilot <- sim_event_times(
  formula = sim_formula,
  data = data_pilot,
  maxt = admin_time
)
lambda_c <- calibrate_censoring(event_times_pilot, censoring_rate, admin_time = admin_time)
cat(sprintf("Censoring distribution lambda in scenario 2: %.4f\n", lambda_c))

trees_2 <- vector("list", length = n_sim)
split_detected_2 <- logical(n_sim)
tau_hat <- numeric(n_sim)
root_p_value_2 <- numeric(n_sim)
observed_censoring_2 <- numeric(n_sim)
errors_2 <- character(0)

seeds <- 959 + seq_len(n_sim)

for (i in seq_len(n_sim)) {
  message("Running simulation ", i, " with seed ", seeds[i])
  set.seed(seeds[i])

  data <- data.frame(
    id = seq_len(n_obs),
    x = rnorm(n_obs),
    z = round(runif(n_obs, min = -1, max = 1), 2)
  )
  sim_data <- sim_surv(
    formula = sim_formula,
    data = data,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
  message("Observed censoring rate: ", mean(sim_data$status == 0))
  observed_censoring_2[i] <- mean(sim_data$status == 0)

  fit <- tryCatch(
    pemtree(
      formula = fit_formula,
      data = sim_data,
      cut = fit_cut,
      maxdepth = 2,
      alpha = 0.05
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    errors_2 <- c(errors_2, conditionMessage(fit))
    trees_2[[i]] <- NA
    split_detected_2[i] <- NA
    root_p_value_2[i] <- NA
    tau_hat[i] <- NA
    next
  }

  tree <- fit$tree
  trees_2[[i]] <- tree
  root <- node_party(tree)
  n_terminal <- length(nodeids(tree, terminal = TRUE))
  split_detected_2[i] <- n_terminal > 1
  tau_hat[i] <- if (n_terminal > 1) {
    breaks_split(split_node(root))
  } else {
    NA_real_
  }
  root_p_value_2[i] <- info_node(root)$p.value
}

power_2 <- mean(split_detected_2, na.rm = TRUE)
cat(sprintf("Power in scenario 2: %.2f", power_2))

# Inspect distribution of tau_hat
hist(tau_hat, breaks = seq(-1, 1, by = 0.05))
abline(v = 0, col = "red")
mean(abs(tau_hat), na.rm = TRUE)
mean(abs(tau_hat) <= 0.1, na.rm = TRUE)
