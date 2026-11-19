## Fluctuation process on PED data vs. theoretical Brownian Bridge

## Settings
library("survival")
library("pammtools")
set.seed(1)

n <- 500
lambda <- 0.1

## Simulate data under H0
sim_ped <- function(n, cut = NULL) {
  T <- rexp(n, lambda)
  C <- runif(n, 0, 20)
  df <- data.frame(id = seq_len(n), time = pmin(T, C), status = as.integer(T <= C))
  ped_args <- list(formula = Surv(time, status) ~ 1, data = df)
  if (!is.null(cut)) ped_args$cut <- cut
  ped <- do.call(as_ped, ped_args)
  ped[order(ped$tend), ]
}

score_process <- function(ped) {
  t_ij <- exp(ped$offset)
  lam <- sum(ped$ped_status) / sum(t_ij)
  psi <- t_ij * lam - ped$ped_status
  W <- cumsum(psi) / sqrt(sum(psi^2))
  end <- which(c(diff(ped$tend) != 0, TRUE))
  list(
    standard = data.frame(
      t = seq_along(W) / length(W),
      W = W
    ),
    time_at_risk = data.frame(
      t = cumsum(t_ij)[end] / sum(t_ij),
      W = W[end]
    )
  )
}

bridge <- function(m = 1000) {
  t <- seq_len(m) / m
  B <- cumsum(rnorm(m, sd = sqrt(1 / m)))
  data.frame(t = t, W = B - t * B[m])
}

## Plot 1: One data set next to theoretical bridges
sp <- score_process(sim_ped(n = 500))
par(mfrow = c(1, 2))
for (v in c("standard", "time_at_risk")) {
  plot(NA, xlim = c(0, 1), ylim = c(-2, 2), xlab = "t", ylab = "W(t)", main = v)
  for (b in 1:20) lines(bridge(), col = "grey80")
  curve(1.96 * sqrt(x * (1 - x)), add = TRUE, lty = 2)
  curve(-1.96 * sqrt(x * (1 - x)), add = TRUE, lty = 2)
  lines(sp[[v]], col = c(standard = "red", time_at_risk = "blue")[v], lwd = 2)
}

## Plot 2: Many data sets, variance an fixed points t
## For convergence against brownian bridge Var(W(t)) = t(1 - t) must hold
R <- 500
t_grid <- seq(0.05, 0.95, by = 0.05)

values <- replicate(R, {
  sp <- score_process(sim_ped(n = 500))
  values_at <- function(d) d$W[findInterval(t_grid, d$t)]
  c(values_at(sp$standard), values_at(sp$time_at_risk))
})
var_standard <- apply(values[seq_along(t_grid), ], 1, var)
var_risk_time <- apply(values[-seq_along(t_grid), ], 1, var)

par(mfrow = c(1, 1))
curve(x * (1 - x), from = 0, to = 1, lwd = 2, ylim = c(0, 0.3), xlab = "t", ylab = "Var(W(t))")
points(t_grid, var_standard, pch = 19, col = "red")
points(t_grid, var_risk_time, pch = 19, col = "blue")
legend(
  "topright",
  bty = "n",
  pch = c(NA, 19, 19),
  lwd = c(2, NA, NA),
  col = c("black", "red", "blue"),
  legend = c("Brownian Bridge", "Standard", "Time-at-risk")
)
