### setup
# packages
source("setup.R")
# simulation functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
# pemtree_splines function
source("R/pemtree_splines.R")
# methods for object of class pemtree
source("R/pemtree.R")

#### Helper functions #########################################################

# Function for plotting baseline hazard
plot_basehaz <- function(hazard_fun) {
  t_vals <- seq(0, 10, by = 0.01)
  plot(t_vals, exp(hazard_fun(t_vals)),
    type = "l",
    xlab = "t", ylab = expression(lambda(t))
  )
}


## DGP 1 #######################################################################
## Two baseline functions with covariates ######################################
################################################################################

# log(lambda(t | x_i)) = beta_0 + f_0(t) * I[t < tau] + g_0(t) * I[t >= tau] +
# beta_1 * x_1 + beta_2 * x_2

# Define baseline
log_hazard <- function(t, tau = 4) {
  ifelse(
    t < tau,
    -2 + 0.5 * t - 0.15 * t^2,
    -1.6 + 0.5 * (sin(t) - sin(4))
  )
}

# Plot baseline hazard
plot_basehaz(log_hazard)

# Simulate data
n <- 2000
covariate_spec <- list(
  x1 = list(
    distfun = rnorm,
    mean = 1,
    sd = 3
  ),
  x2 = list(
    distfun = rbinom,
    prob = 0.5,
    size = 1
  )
)
formula <- ~ log_hazard(t) + 0.08 * x1 - 0.1 * x2
censoring_rate <- 0.2
admin_time <- 10L
cut <- seq(0, admin_time, by = 0.1)

set.seed(1746)
sim_data <- sim_pam(
  n = n,
  formula = formula,
  covariate_spec = covariate_spec,
  censoring_rate = censoring_rate,
  admin_time = admin_time,
  cut = cut
)
prop.table(table(sim_data$censoring_reason))

# Fit a pemtree with splines
tree_fit_1 <- pemtree_splines(
  formula = Surv(time, status) ~ bs(tend, df = 4, degree = 3) + x1 + x2 +
    offset(offset) | tend,
  data = sim_data,
  cut = cut,
  min_events = 20,
  maxdepth = 2,
  alpha = 0.05
)

tree_fit_1

# Create node specific baseline hazard predictions
tree1 <- tree_fit_1$tree
terminal_ids <- nodeids(tree1, terminal = TRUE)
node_models <- nodeapply(
  tree1,
  ids = terminal_ids,
  FUN = function(node) info_node(node)$object
)

newdata <- data.frame(
  tend = seq(min(sim_data$time), max(sim_data$time), length.out = 500),
  offset = 0, x1 = 0, x2 = 0
)
newdata$node_id <- predict(tree1, newdata = newdata, type = "node")
newdata$basehaz_hat <- NA_real_
for (j in seq_along(terminal_ids)) {
  rows <- terminal_ids[j] == newdata$node_id
  newdata$basehaz_hat[rows] <- suppressWarnings(predict(
    node_models[[j]],
    newdata = newdata[rows, ],
    type = "response"
  ))
}
newdata$basehaz_true <- exp(log_hazard(newdata$tend))

# Plot estimate vs true baseline hazard
ggplot(newdata, aes(x = tend)) +
  geom_line(aes(y = basehaz_hat, linetype = "estimate")) +
  geom_line(aes(y = basehaz_true, linetype = "true")) +
  labs(
    x = "time",
    y = expression(lambda(t)),
    linetype = NULL
  ) +
  theme_bw()

## DGP 2 #######################################################################
## Time-varying covariate effects ##############################################
################################################################################

# log(lambda(t | x)) = f_0(t) + beta_1 * x_1 + beta_2 * x_1 * I[t >= tau]
log_hazard <- function(t) {
  -2.5 + 1.4 * dgamma(t, shape = 8, rate = 2)
}

# Plot baseline hazard
plot_basehaz(log_hazard)

# Simulate data
n <- 2000
covariate_spec <- list(x1 = list(distfun = rnorm, mean = 1, sd = 2))
tau <- 5
formula <- ~ log_hazard(t) + 0.3 * x1 + 0.5 * x1 * I(t >= 5)
censoring_rate <- 0.2
admin_time <- 10
cut <- seq(0, admin_time, by = 0.01)

set.seed(1912)
sim_data <- sim_pam(
  n = n,
  formula = formula,
  covariate_spec = covariate_spec,
  censoring_rate = censoring_rate,
  admin_time = admin_time,
  cut = cut
)

# Fit a pemtree with splines
tree_fit_2 <- pemtree_splines(
  formula = Surv(time, status) ~ bs(tend, df = 5, degree = 3) + x1 +
    offset(offset) | tend,
  data = sim_data,
  cut = seq(0, max(sim_data$time), by = 0.1),
  min_events = 30,
  maxdepth = 2,
  alpha = 0.05
)

tree_fit_2
tree2 <- tree_fit_2$tree

# Extract time split
tau_hat <- breaks_split(split_node(node_party(tree2)))
tau_hat
# Extract coefficients for x1
coef_tree2 <- coef(tree2, drop = FALSE)
coef_tree2_x1 <- coef_tree2[, "x1"]
coef_tree2_x1
