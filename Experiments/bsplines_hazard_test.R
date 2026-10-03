### setup
source("setup.R")
## Load required functions
source("R/pemtree_splines.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")

## Data-generating process 1
log_hazard <- function(t, tau = 4, delta) {
  smooth_baseline <- function(t) {
    -3 + 0.5 * sin(pi * t / 4) + 0.05 * t
  }

  smooth_baseline(t) + delta * (t >= tau)
}

# Plot baseline hazard
t_vals <- seq(0, 8, by = 0.01)
plot(t_vals, exp(log_hazard(t_vals, delta = 1)),
  type = "l",
  xlab = "t", ylab = expression(lambda(t))
)

# Simulate data
n <- 2000
formula <- ~ log_hazard(t, delta = 1)
censoring_rate <- 0.3
cut <- seq(0, 8, by = 0.1)
set.seed(2103)
sim_data <- sim_pam(
  n = n,
  formula = formula,
  covariate_spec = NULL,
  censoring_rate = censoring_rate,
  admin_time = max(cut),
  cut = cut
)
# Inspect censoring rate
prop.table(table(sim_data$status))
prop.table(table(sim_data$censoring_reason))[c("random", "administrative")]

# Fit the pemtree
tree_fit <- pemtree_splines(
  formula = Surv(time, status) ~ bs(tend, degree = 3, df = 5) +
    offset(offset) | tend,
  data = sim_data,
  cut = cut,
  min_events = 30,
  maxdepth = 2,
  alpha = 0.05
)

tree_fit$tree

## Reconstruct the node-specific hazard estimate
tree <- tree_fit$tree
terminal_ids <- nodeids(tree, terminal = TRUE)
node_models <- nodeapply(
  tree,
  ids = terminal_ids,
  FUN = function(node) info_node(node)$object
)
newdata <- data.frame(
  tend = seq(min(sim_data$time), max(sim_data$time), length.out = 500),
  offset = 0
)
newdata$node_id <- predict(tree, newdata = newdata, type = "node")
newdata$hazard_hat <- NA_real_
for (j in seq_along(terminal_ids)) {
  rows <- newdata$node_id == terminal_ids[[j]]
  newdata$hazard_hat[rows] <- predict(
    tree,
    newdata = newdata[rows, ],
    type = "response"
  )
}
newdata$hazard_true <- exp(log_hazard(newdata$tend, delta = 1))

ggplot(newdata, aes(x = tend)) +
  geom_line(aes(y = hazard_hat, linetype = "Estimate")) +
  geom_line(aes(y = hazard_true, linetype = "True")) +
  labs(x = "t", y = expression(lambda(t))) +
  theme_bw()

## Data-generating process 2
log_hazard <- function(t, tau = 6) {
  ifelse(
    t < tau,
    -3 + 2 * dgamma(t, shape = 8, rate = 2),
    -1.5
  )
}
# Plot baseline hazard
t_vals <- seq(0, 10, by = 0.01)
plot(t_vals, exp(log_hazard(t_vals, tau = 6)),
  type = "l",
  xlab = "t", ylab = expression(lambda(t))
)

covariate_spec <- list(x1 = list(distfun = runif, min = 0, max = 4))
formula <- ~ log_hazard(t, tau = 6) + sqrt(x1)
cut <- seq(0, 10, by = 0.1)
n <- 2000

set.seed(1111)
sim_data <- sim_pam(
  n = n,
  formula = formula,
  covariate_spec = covariate_spec,
  censoring_rate = 0.2,
  admin_time = 10,
  cut = cut
)

# Fit the tree
tree_fit <- pemtree_splines(
  formula = Surv(time, status) ~ bs(tend, df = 5, degree = 3) + bs(x1, df = 5, degree = 3) +
    offset(offset) | tend,
  data = sim_data,
  cut = cut,
  min_events = 30,
  maxdepth = 2,
  alpha = 0.05
)

tree_fit$tree
