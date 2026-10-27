############################################################
## Analysis of pamtree #####################################

### setup
source("setup.R")
source("R/sim_surv.R")

## Simulate data without a level shift
shape <- 1.5
rate <- 0.05
log_hazard <- function(t) log(rate) + log(shape) + (shape - 1) * log(t)

formula <- ~ log_hazard(t)

admin_time <- 10
censoring_rate <- 0.3

# Calibrate lambda for censoring distribution
set.seed(12)
n_pilot <- 10000
data_pilot <- data.frame(id = seq_len(n_pilot))
event_times_pilot <- sim_event_times(formula, data_pilot)
lambda_c <- calibrate_censoring(event_times_pilot, censoring_rate, admin_time)

# Simulate time to event data
set.seed(22)
n <- 1000
data <- data.frame(id = seq_len(n))
sim_data <- sim_surv(formula, data, admin_time = admin_time, lambda_c = lambda_c)
mean(sim_data$status == 0)

# Transform to PED data
ped <- as_ped(Surv(time, status) ~ ., data = sim_data)

# Fit a gam with all observations
fit_gam <- gam(
  ped_status ~ s(tend, bs = "ps") + offset(offset),
  family = poisson(link = "log"),
  data = ped,
  method = "REML"
)

# Extract scores like pamtree does
X <- predict(fit_gam, type = "lpmatrix")
raw_scores <- sweep(X, 1, ped$ped_status - fitted(fit_gam), FUN = "*")
scores <- sweep(raw_scores, 2, colMeans(raw_scores), FUN = "-")
