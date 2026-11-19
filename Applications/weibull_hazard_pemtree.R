source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_exact.R")
### Weibull Hazard on exposure-scaled pemtree -----------------------------
# Investigate how well a split can be found with a Weibull baseline hazard

# Function to generate data
generate_data <- function(n, delta, tau, lambda_c) {
  sim_surv(
    ~ log(0.05) + log(1.5) + 0.5 * log(t) + delta * (t >= tau),
    data = data.frame(id = seq_len(n)),
    breaks = tau,
    admin_time = 10,
    lambda_c = lambda_c
  )
}

# Oracle model
set.seed(1)
data_oracle <- generate_data(n = 1000, delta = 0.5, tau = 3, lambda_c = 0.1)
ped <- as_ped(data_oracle, Surv(time, status) ~ 1)
fit_orcacle <- glm(
  ped_status ~ log(tend) + I(tend >= 3),
  family = poisson(),
  offset = offset,
  data = ped
)
summary(fit_orcacle)

# Predict with new data
t_grid <- seq(0.1, 10, length.out = 500)
newdata <- data.frame(tend = t_grid, offset = 0)
newdata$hazard_hat <- predict(fit_orcacle, newdata = newdata, type = "response")
newdata$hazard_true <- exp(log(0.05) + log(1.5) + 0.5 * log(newdata$tend) +
  0.5 * (newdata$tend >= 3))

ggplot(newdata, aes(x = t_grid)) +
  geom_line(aes(y = hazard_hat, linetype = "Estimate")) +
  geom_line(aes(y = hazard_true, linetype = "True")) +
  labs(x = "t", y = expression(lambda(t)), linetype = NULL) +
  theme_bw()

# Fit exposure-scaled pemtree once
# Simulate data
set.seed(2)
data <- generate_data(n = 500, delta = 0.5, tau = 3, lambda_c = 0.1)
# Check censoring
mean(data$status == 0)

# Fit
ped <- as_ped(data, Surv(time, status) ~ 1)
fit <- pem_tree(
  ped = ped,
  alpha = 0.05,
  pi0 = 0.1,
  R = 1e4,
  minevents = 20,
  maxdepth = 2
)
fit

# Repeat procedure 1000 times
