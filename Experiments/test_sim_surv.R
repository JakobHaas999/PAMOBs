####################################################
## Example DGPs for sim_surv #######################
####################################################

source("setup.R")
source("R/sim_surv.R")

## Scenario 1: Constant hazard ------------------------------------------------------
log_hazard <- ~ -2
data <- data.frame(id = seq_len(1000))

set.seed(1)
sim_data_1 <- sim_surv(
  formula = log_hazard, data = data, censoring_rate = 0.3
)
mean(sim_data_1$status == 0)

ped <- pammtools::as_ped(
  Surv(time, status) ~ 1,
  data = sim_data_1, cut = seq(0, ceiling(max(sim_data_1$time)), by = 1)
)
data.table::setDT(ped)

hazard_est <- ped[, .(
  events = sum(ped_status), exposure = sum(exp(offset)), n_risk = .N
), by = .(tstart, tend)]
hazard_est[, hazard := events / exposure]
head(hazard_est)


## Scenario 2: Weibull Hazard without structural change ------------------------------
shape <- 1.5
rate <- 0.05

log_hazard <- ~ log(rate) + log(shape) + (shape - 1) * log(t)
data <- data.frame(id = seq_len(3000))

set.seed(2)
sim_data_2 <- sim_surv(
  formula = log_hazard, data = data, censoring_rate = 0.3
)

ped_2 <- pammtools::as_ped(
  Surv(time, status) ~ 1,
  data = sim_data_2, cut = seq(0, ceiling(max(sim_data_2$time)), by = 1)
)
data.table::setDT(ped_2)

hazard_est_2 <- ped_2[, .(
  events = sum(ped_status), exposure = sum(exp(offset)), n_risk = .N
), by = .(tstart, tend)]
hazard_est_2[, hazard := events / exposure]
head(hazard_est_2)

truth <- data.frame(t = seq(0, 10, length.out = 500))
truth$hazard <- rate * shape * truth$t^(shape - 1)

ggplot() +
  geom_segment(
    data = hazard_est_2[tstart < 10],
    aes(x = tstart, xend = tend, y = hazard, yend = hazard),
    colour = "steelblue", linewidth = 0.9
  ) +
  geom_line(data = truth, aes(x = t, y = hazard)) +
  labs(x = "time", y = expression(lambda(t))) +
  theme_bw()

## ## Scenario 2: Weibull Hazard with structural change ---------------------------------
shape <- 1.5
rate <- 0.05
delta <- log(3)
tau <- 5

log_hazard <- ~ log(rate) + log(shape) + (shape - 1) * log(t) + delta * (t >= tau)
data <- data.frame(id = seq_len(3000))

set.seed(3)
sim_data_3 <- sim_surv(
  formula = log_hazard, data = data, breaks = tau, censoring_rate = 0.3
)

truth <- data.frame(t = sort(unique(c(seq(0, 10, by = 0.01), tau))))
truth$hazard <- with(truth, rate * shape * t^(shape - 1) * exp(delta))

cuts <- sort(unique(c(
  seq(0, ceiling(max(sim_data_3$time)), by = 1), tau
)))
ped_3 <- pammtools::as_ped(
  formula = Surv(time, status) ~ 1, data = sim_data_3, cut = cuts
)
data.table::setDT(ped_3)
hazard_est_3 <- ped_3[, .(
  events = sum(ped_status), exposure = sum(exp(offset)), n_risk = .N
), by = .(tstart, tend)]
hazard_est_3[, hazard := events / exposure]

ggplot() +
  geom_segment(
    data = hazard_est_3[tstart < 10],
    aes(x = tstart, xend = tend, y = hazard, yend = hazard),
    colour = "steelblue", linewidth = 0.9
  ) +
  geom_line(data = truth, aes(x = t, y = hazard)) +
  labs(x = "time", y = expression(lambda(t))) +
  theme_bw()
