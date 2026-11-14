############################################################################
## Complex tree DGP. different baseline between subgroups and level shift ##

## setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

## DGP:
#       root
#     /     \
# z = 0     z = 1
#           /   \
#         t <= 3 t > 3

beta0_0 <- -2
beta0_1 <- -1.5

tau <- 3
delta <- 0.5

sim_formula <- ~ beta0_0 * (z == 0) + beta0_1 * (z == 1) + delta * (t > tau) * (z == 1)

data <- sim_surv(
  sim_formula, data.frame(id = seq_len(1000), z = rbinom(1000, size = 1, prob = 0.5)),
  breaks = tau, admin_time = 10, lambda_c = 0.1
)
mean(data$status == 0)

ped <- as_ped(
  data, Surv(time, status) ~ z,
  cut = c(0, tau, 10)
)
ped$z_late <- as.integer(ped$z == 1 & ped$tend > tau)

# Fit oracle model
oracle <- glm(
  ped_status ~ 0 + factor(z) + z_late,
  offset = offset,
  family = poisson(link = "log"),
  data = ped
)
summary(oracle)

nd <- expand.grid(tend = seq(0.1, 10, length.out = 500), z = factor(c(0, 1)))
nd$z_late <- as.integer(nd$z == 1 & nd$tend > tau)
nd$offset <- 0
nd$hazard_hat <- predict(oracle, newdata = nd, type = "response")
nd$hazard_true <- exp(ifelse(nd$z == 0, beta0_0, beta0_1 + delta * (nd$tend > tau)))

ggplot(nd, aes(x = tend, colour = z)) +
  geom_line(aes(y = hazard_hat)) +
  geom_line(aes(y = hazard_true), linetype = "dashed") +
  labs(x = "t", y = expression(lambda(t))) +
  theme_bw()

# Fit pemtree
tree_fit <- pemtree(
  Surv(time, status) ~ offset(offset) | tend + z,
  data = data,
  cut = seq(0, 10, by = 0.1),
  min_events = 20,
  alpha = 0.05
)

plot(tree_fit$tree)

nd$hazard_hat_tree <- predict(tree_fit$tree, newdata = nd, type = "response")
ggplot(nd, aes(x = tend, colour = z)) +
  geom_line(aes(y = hazard_hat_tree)) +
  geom_line(aes(y = hazard_true), linetype = "dashed") +
  labs(x = "t", y = expression(lambda(t))) +
  theme_bw()
