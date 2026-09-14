# NOTE:
# PED rows from the same subject are not independent.
# This analysis is currently exploratory and is used to assess whether
# a temporal structural break is visible to MOB at all.

# ==============================================================================
# Detecting Temporal Changes in the Baseline Hazard with MOB
#
# The aim of this script is to investigate whether model-based recursive
# partitioning (MOB) can detect a structural change in the baseline hazard.
#
# Survival times are simulated from a piecewise constant hazard of the form
#
#   lambda(t) = lambda_1,  for t <= tau
#             = lambda_2,  for t >  tau,
#
# where tau denotes the true temporal change point.
#
# The simulated survival data are transformed into piecewise exponential data
# (PED), and MOB is applied with time as the partitioning variable. The goal is
# to assess whether MOB detects the change in the hazard and recovers the true
# change point tau.
# ==============================================================================

## load required files
source("setup.R")
source("R/simulation_pem.R")
source("R/simulation_censoring.R")

### Simulate data with one structural change point
set.seed(12)
n <- 10000
lambdas <- c(0.2, 0.7)
cuts <- 4

time <- sim_pem(n = n, lambdas = lambdas, cuts = cuts)

#### No censoring ####################################################
######################################################################
surv_data <- data.frame(time = time, status = 1)

# Transform data to PED
ped <- pammtools::as_ped(
  formula = Surv(time, status) ~ 1,
  data = surv_data,
  cut = seq(0, max(surv_data$time), by = 0.25)
)

## Compare normal PEM with MOB
# 1) Fit PEM
pem_fit <- glm(
  formula = ped_status ~ interval,
  family = poisson(link = "log"),
  data = ped,
  offset = offset
)
newdata_pem <- data.frame(offset = 0, interval = levels(ped$interval))
hazard_pem <- predict(
  pem_fit,
  newdata = newdata_pem,
  type = "response"
)

# 2) Fit MOB
tree_fit <- partykit::glmtree(
  formula = ped_status ~ 1 | tend,
  family = poisson(link = "log"),
  data = ped,
  offset = ped$offset,
  alpha = 0.05
)
newdata_tree <- data.frame(
  tend = seq(min(ped$tend), max(ped$tend), by = 0.25),
  offset = 0
)
nodes <- predict(
  tree_fit,
  newdata = newdata_tree,
  type = "node"
)
terminal_ids <- partykit::nodeids(
  tree_fit,
  terminal = TRUE
)
node_models <- partykit::nodeapply(
  tree_fit,
  ids = terminal_ids,
  FUN = function(node) node$info$object
)
node_hazards <- vapply(node_models, function(mod) {
  exp(coef(mod))["(Intercept)"]
}, numeric(1))

hazard_tree <- node_hazards[match(nodes, terminal_ids)]

## Plot differences
hazard_df <- cbind.data.frame(
  tend = newdata_tree$tend,
  hazard_pem,
  hazard_tree
)
hazard_df$hazard_true <- ifelse(hazard_df$tend <= 4, lambdas[[1]], lambdas[[2]])

ggplot(hazard_df, aes(x = tend)) +
  geom_step(aes(y = hazard_pem, colour = "pem estimate")) +
  geom_step(aes(y = hazard_tree, colour = "tree estimate")) +
  geom_step(aes(y = hazard_true, colour = "true")) +
  geom_vline(xintercept = cuts, linetype = "dashed") +
  labs(
    x = "Time",
    y = "Hazard",
    linetype = NULL,
    title = "True and Estimates for Baseline Hazard"
  ) +
  theme_bw()


#### With censoring ##################################################
######################################################################

# Simulate censoring times
set.seed(13)
censoring_rate <- 0.6
cens_surv_data <- sim_censoring(times = time, censoring_rate = censoring_rate)

# Transform data
ped_cens <- pammtools::as_ped(
  formula = Surv(time, status) ~ 1,
  data = cens_surv_data,
  cut = seq(0, max(cens_surv_data$time), by = 0.25)
)

## Fit a MOB
tree_cens <- partykit::glmtree(
  formula = ped_status ~ 1 | tend,
  family = poisson(link = "log"),
  data = ped_cens,
  offset = ped_cens$offset,
  alpha = 0.05
)

# Plot estimated hazard vs. true hazard
nodes_cens <- predict(
  tree_cens,
  newdata = data.frame(
    offset = 0,
    tend = seq(min(ped_cens$tend), max(ped_cens$tend), by = 0.25)
  ),
  type = "node"
)
terminal_ids_cens <- partykit::nodeids(
  tree_cens,
  terminal = TRUE
)
node_models_cens <- nodeapply(
  tree_cens,
  ids = terminal_ids_cens,
  FUN = function(node) node$info$object
)
node_hazards_cens <- vapply(node_models_cens, function(mod) {
  exp(coef(mod))["Intercept"]
}, numeric(1))

hazard_cens <- node_hazards_cens[match(nodes_cens, terminal_ids_cens)]

hazard_df_cens <- cbind.data.frame(
  tend = seq(min(ped_cens$tend), max(ped_cens$tend), by = 0.25),
  hazard_tree = hazard_cens
)
hazard_df_cens$hazard_true <- ifelse(hazard_df_cens$tend < 4, lambdas[[1]], lambdas[[2]])

ggplot(hazard_df_cens, aes(x = tend)) +
  geom_step(aes(y = hazard_tree, colour = "estimate")) +
  geom_step(aes(y = hazard_true, colour = "true")) +
  geom_vline(xintercept = cuts, linetype = "dashed") +
  labs(
    x = "Time",
    y = "Hazard",
    linetype = NULL,
    title = "True and MOB Estimate for Baseline Hazard"
  ) +
  theme_bw()
