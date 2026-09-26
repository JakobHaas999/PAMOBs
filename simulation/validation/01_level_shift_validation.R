source("simulation/scenarios/01_level_shift.R")
tau <- 4
beta0 <- -2
delta <- 1

log_hazard <- function(t) {
  beta0 + delta * (t > tau)
}
t_grid <- seq(0, 8, by = 0.01)

## Simulate data
set.seed(123)
dat <- scenario_level_shift(
  n = 10000,
  delta = delta,
  tau = tau,
  beta0 = beta0
)

ped <- as_ped(
  formula = Surv(time, status) ~ 1,
  data = dat,
  cut = seq(0, 8, by = 0.1)
)

fit <- glm(
  ped_status ~ interval + offset(offset),
  family = poisson(link = "log"),
  data = ped
)


## Compare prediction with truth
# Prediction data for all intervalls
newdata <- data.frame(interval = levels(ped$interval), offset = 0)
hazard_hat <- predict(fit, newdata = newdata, type = "response")
interval_info <- unique(ped[, c("interval", "tend")])
interval_info <- interval_info[
  match(levels(ped$interval), interval_info$interval),
]
interval_info$hazard_hat <- hazard_hat
interval_info$hazard_true <- exp(log_hazard(interval_info$tend))
head(interval_info)

## Plot prediction and truth
ggplot(interval_info, aes(x = tend)) +
  geom_step(aes(y = hazard_hat, col = "estimate")) +
  geom_step(aes(y = hazard_true, col = "truth")) +
  ylab(expression(lambda(t))) +
  xlab("t") +
  theme_bw()
