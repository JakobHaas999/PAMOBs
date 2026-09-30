### setup
source("setup.R")
## Load required functions
source("R/pemtree.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")

## Data-generating process

# Before tau, the log-baseline hazard follows a quadratic function.
# At tau, both its functional form and level change; afterwards,
# the log-baseline hazard is constant at -2.
log_hazard <- function(t) {
  ifelse(
    t < 4,
    -3 + 0.8 * t - 0.08 * t^2,
    -2
  )
}
n <- 1000
formula <- ~ log_hazard(t)
censoring_rate <- 0.4
cut <- seq(0, 8, by = 0.01)

set.seed(12)
sim_data <- sim_pam(
  n = n,
  formula = formula,
  covariate_spec = NULL,
  censoring_rate = censoring_rate,
  cut = cut
)

## Fit the temporal PEM tree
fit_cut <- unique(c(seq(0, max(sim_data$time), by = 0.1), max(sim_data$time)))
tree_fit <- pemtree(
  formula = Surv(time, status) ~ poly(tend, 2, raw = TRUE) + offset(offset) |
    tend,
  data = sim_data,
  cut = fit_cut,
  min_events = 30,
  maxdepth = 2,
  alpha = 0.05
)
# Each node models a quadratic log-hazard. MOB searches for a time point
# at which the intercept, slope, or curvature becomes unstable.
# The offset converts expected interval counts into a hazard model.
tree_fit

## Reconstruct the node-specific hazard estimates

# Assign each grid point to its terminal node and extract the
# corresponding polynomial coefficients.

# For each time point, compute
# eta(t) = beta0 + beta1 * t + beta2 * t^2,
# followed by lambda(t) = exp(eta(t)).
# The exposure offset is omitted because the target is the hazard itself.
newdata <- data.frame(tend = seq(min(sim_data$time), max(sim_data$time), length.out = 500))
newdata$node_id <- predict(tree_fit$tree, newdata = newdata, type = "node")
terminal_nodes <- nodeids(tree_fit$tree, terminal = TRUE)
coef_nodes <- coef(tree_fit$tree, node = terminal_nodes, drop = FALSE)
rownames(coef_nodes) <- terminal_nodes

x_grid <- cbind(1, newdata$tend, newdata$tend^2)
beta_grid <- coef_nodes[as.character(newdata$node_id), , drop = FALSE]
newdata$hazard_hat <- exp(rowSums(x_grid * beta_grid))
newdata$hazard_true <- exp(log_hazard(newdata$tend))


## Plot hazard estimate vs. true hazard
ggplot(newdata, aes(x = tend)) +
  geom_line(aes(y = hazard_hat, linetype = "Estimate")) +
  geom_line(aes(y = hazard_true, linetype = "True")) +
  labs(x = expression(t), y = expression(lambda(t)), linetype = NULL) +
  theme_bw()
