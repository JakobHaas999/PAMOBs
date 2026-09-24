source("setup.R")
source("R/simulation_pam.R")

######## Experiment:
### We want to find a structural break in the baseline hazard
### using a fixed spline basis, no penalization
### We want to find the break through model based partitioning

## Simulate some data
set.seed(12)
log_baseline <- function(t) {
  -2 + 0.15 * t + ifelse(t > 4, 1, 0)
}
formula <- ~ log_baseline(t)
cut <- seq(0, 8, by = 0.25)
surv_data <- sim_pam(
  n = 10000,
  formula = formula,
  covariate_spec = NULL,
  baseline = log_baseline,
  cut = cut
)

ped <- as_ped(Surv(time, status) ~ 1, data = surv_data, cut = cut)

fit <- gam(
  formula = ped_status ~ s(tend, bs = "ps"),
  family = poisson(link = "log"),
  data = ped,
  offset = ped$offset,
  method = "REML"
)

X <- predict(fit, type = "lpmatrix")
mu <- fitted(fit)
scores <- X * (ped$ped_status - mu)
scores_c <- scale(scores, center = TRUE, scale = FALSE)

scores_c <- scores_c[ord, , drop = FALSE]
scores_by_time <- rowsum(
  scores_c,
  group = ped$tend,
  reorder = TRUE
)
time <- as.numeric(rownames(scores_by_time))

scores_by_time[nrow(scores_by_time), ]

matplot(time, scores_by_time,
  type = "l", lty = 1,
  xlab = "tend", ylab = "Cumulative centered scores"
)
abline(h = 0, lty = 2)
abline(v = 4, lty = 2)
