source("setup.R")
source("simulation/scenarios/02_level_shift_splines.R")

## Global parameters
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")

basehaz <- function(t) {
  -2.5 + 0.25 * t - 0.025 * t^2
}
delta <- 0.8
tau <- 4
censoring_rate <- 0.2
admin_time <- 10

## Without covariates -----------------------------------------------------------

## Simulate data
set.seed(123)
dat <- scenario_level_shift_splines(
  n = 1000,
  delta = delta,
  tau = tau,
  basehaz_fun = basehaz,
  censoring_rate = censoring_rate,
  admin_time = admin_time
)

ped <- as_ped(
  formula = Surv(time, status) ~ 1,
  data = dat,
  seq(0, admin_time, by = 0.1)
)

fit <- glm(
  formula = ped_status ~ bs(tend, df = 5, degree = 3) + offset(offset),
  family = poisson(link = "log"),
  data = ped
)

# Compare prediction with true hazard
newdata <- data.frame(tend = seq(min(ped$tend), max(ped$tend), length.out = 500), offset = 0)
newdata$hazard_hat <- predict(fit, newdata = newdata, type = "response")
newdata$hazard_true <- exp(basehaz(newdata$tend) + delta * (newdata$tend >= tau))

ggplot(newdata, aes(x = tend)) +
  geom_line(aes(y = hazard_hat, linetype = "Estimate")) +
  geom_line(aes(y = hazard_true, linetype = "True")) +
  labs(x = "t", y = expression(lambda(t)), linetype = NULL) +
  theme_bw()

## With covariates -----------------------------------------------------------

covariate_spec <- list(
  x1 = list(
    distfun = rnorm
  ),
  x2 = list(
    distfun = runif,
    min = -2,
    max = 2
  ),
  x3 = list(
    distfun = rbinom,
    prob = 0.5,
    size = 1
  )
)
covariate_effects <- c(x1 = 0.7, x2 = 1.2, x3 = 0.08)

# Simulate data
set.seed(234)
dat <- scenario_level_shift_splines(
  n = 10000,
  delta = delta,
  tau = tau,
  basehaz_fun = basehaz,
  censoring_rate = censoring_rate,
  covariate_spec = covariate_spec,
  covariate_effects = covariate_effects,
  admin_time = admin_time
)

ped <- as_ped(
  Surv(time, status) ~ x1 + x2 + x3,
  data = dat,
  cut = seq(0, admin_time, by = 0.1)
)

fit <- glm(
  formula = ped_status ~ bs(tend, df = 5, degree = 3) + x1 + x2 + x3 + offset(offset),
  family = poisson(link = "log"),
  data = ped
)

# Compute prediction with true hazard
newdata <- data.frame(
  tend = seq(min(ped$tend), max(ped$tend), length.out = 500),
  offset = 0, x1 = 0, x2 = 0, x3 = 0
)
newdata$hazard_hat <- predict(fit, newdata = newdata, type = "response")
newdata$hazard_true <- exp(basehaz(newdata$tend) + delta * (newdata$tend >= tau))

ggplot(newdata, aes(x = tend)) +
  geom_line(aes(y = hazard_hat, linetype = "Estimate")) +
  geom_line(aes(y = hazard_true, linetype = "True")) +
  labs(x = "t", y = expression(lambda(t)), linetype = NULL) +
  theme_bw()
