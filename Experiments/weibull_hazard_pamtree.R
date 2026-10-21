###################################################################
# Pamtree for Weibull distributed survival times ------------------
###################################################################

source("setup.R")
source("R/sim_surv.R")
source("R/pamtree.R")

## Helper functions -----------------------------------------------
predict_hazard <- function(object, newdata, ...) {
  newdata$offset <- 0
  as.numeric(
    predict(object$tree, newdata = newdata, type = "response")
  )
}


## Analysis -------------------------------------------------------
## 3 scenarios
# 1) Null scenario: Hazard without level shift -> No split.
# 2) Time-varying covariate effect. The effect of a covariate changes at some time point tau.
#    Can we recover a split point here?
# 3) Combine a weibull Hazard with a exponential hazard with a level shift as well
#    Can we recover a split here as well?

## Scenario 1 -----------------------------------------------------

# Parameters
shape <- 1.5
rate <- 0.05

formula <- ~ log(rate) + log(shape) + (shape - 1) * log(t)
data <- data.frame(id = 1:1000)

admin_time <- 10
censoring_rate <- 0.3
n_sim <- 20

split_detected <- logical(n_sim)
tau_hats <- numeric(n_sim)
models <- list()
p_vals <- numeric(n_sim)

seeds <- 2021L + seq_len(n_sim)

for (i in seq_len(n_sim)) {
  cat("\nRunning simulation", i, "with seed", seeds[i], "\n")
  set.seed(seeds[i])
  sim_data <- sim_surv(
    formula = formula,
    data = data,
    admin_time = admin_time,
    censoring_rate = censoring_rate
  )
  cat("Observed censoring rate:", mean(sim_data$status == 0), "\n")

  fit_cut <- seq(0, ceiling(max(sim_data$time)), by = 0.05)
  fit_tree <- pamtree(
    formula = Surv(time, status) ~ s(tend, bs = "ps", k = 8) + offset(offset) | tend,
    data = sim_data,
    cut = fit_cut,
    gam_method = "ML",
    sp = 1,
    alpha = 0.05,
    maxdepth = 2
  )

  n_terminal_nodes <- length(nodeids(fit_tree$tree, terminal = TRUE))
  root <- node_party(tree_fit$tree)
  p_vals[i] <- info_node(root)$test["p.value", ]
  tau_hat <- if (n_terminal_nodes == 2) {
    root <- node_party(fit_tree$tree)
    breaks_split(split_node(root))
  } else {
    NA_real_
  }
  tau_hats[i] <- tau_hat
  split_detected[i] <- n_terminal_nodes > 1
  models[[length(models) + 1]] <- fit_tree
}

mean(split_detected)

# Plot estimated hazards
t_grid <- seq(0, admin_time, length.out = 500)
hazard_true <- rate * shape * (t_grid)^(shape - 1)
newdata <- data.frame(tend = t_grid)
hazard_estimates <- lapply(models, function(model) predict_hazard(model, newdata))

plot(
  t_grid, hazard_true,
  type = "l",
  col = "black",
  lwd = 3,
  ylim = range(c(hazard_true, unlist(hazard_estimates))),
  xlab = "t",
  ylab = expression(lambda(t)),
  main = "Estimated hazard"
)
for (i in seq_along(hazard_estimates)) {
  lines(
    t_grid, hazard_estimates[[i]],
    col = adjustcolor("steelblue", alpha.f = 0.3),
    lwd = 1.5
  )
}
