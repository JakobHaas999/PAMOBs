source("setup.R")
source("R/simulation_pam.R")

######## Experiment:
### We want to find a structural break in the baseline hazard
### using a fixed spline basis, no penalization
### We want to find the break through model based partitioning

## Simulate some data
set.seed(12)
log_baseline <- function(t) {
  -2 + sin(t) + ifelse(t > 4, 3, 0)
}
formula <- ~ log_baseline(t)
cut <- seq(0, 8, by = 0.25)
censoring_rate <- 0.4

surv_data <- sim_pam(
  n = 10000,
  formula = formula,
  covariate_spec = NULL,
  baseline = log_baseline,
  censoring_rate = censoring_rate,
  cut = cut
)
head(surv_data)
prop.table(table(surv_data$status))

ped <- as_ped(Surv(time, status) ~ 1, data = surv_data, cut = cut)

## B-Spline Basis
B <- splines::bs(ped$tend, df = 3, degree = 3, intercept = FALSE, Boundary.knots = c(0, 8))
colnames(B) <- paste0("B", 1:3)
ped <- cbind.data.frame(ped, B)

tree <- palmtree::palmtree(
  formula = ped_status ~ 1 | B1 + B2 + B3 + offset(offset) | tend,
  data = ped,
  family = poisson()
)

# Evaluate hazard
t_grid <- seq(min(ped$tend), max(ped$tend), length.out = 500)
B_grid <- predict(B, newx = t_grid)
colnames(B_grid) <- paste0("B", 1:3)
newdata <- data.frame(tend = t_grid, offset = 0, B_grid)
hazard_hat <- predict(
  object = tree,
  newdata = newdata,
  type = "response"
)

plot(t_grid, hazard_hat,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = "Hazard"
)
lines(t_grid, exp(log_baseline(t_grid)),
  lwd = 2,
  lty = 2
)
abline(v = 4, lty = 3)
legend(
  "topright",
  legend = c("PALM tree", "True hazard"),
  lty = c(1, 2),
  lwd = 2
)


# p <- poly(ped$tend, degree = 3, raw = FALSE)
# colnames(p) <- paste0("p", 1:3)
# ped <- cbind.data.frame(ped, p)
#
# tree <- glmtree(
#   formula = ped_status ~ p1 + p2 + p3 + offset(offset) | tend,
#   data = ped,
#   family = poisson(link = "log")
# )
#
# t_grid <- seq(min(ped$tend), max(ped$tend), length.out = 500)
# P_grid <- predict(p, newdata = t_grid)
# colnames(P_grid) <- paste0("p", 1:3)
# newdata <- data.frame(
#   tend = t_grid,
#   p1 = P_grid[, 1],
#   p2 = P_grid[, 2],
#   p3 = P_grid[, 3],
#   offset = 0
# )
#
# hazard_hat <- predict(tree, newdata = newdata, type = "response")
