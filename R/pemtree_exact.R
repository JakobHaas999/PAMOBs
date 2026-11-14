agg_ped <- function(ped) {
  tij <- exp(ped$offset)
  M <- rowsum(
    cbind(T = tij, D = ped$ped_status, n = 1), ped$tend,
    reorder = TRUE
  )
  list(
    tend = sort(unique(ped$tend)),
    D_j = M[, "D"],
    T_j = M[, "T"],
    n_j = M[, "n"],
    tij = tij
  )
}

suplm_time <- function(ped, pi_0 = 0.1, R = 1e4) {
  a <- agg_ped(ped)
  D <- sum(a$D_j)
  Tt <- sum(a$T_j)
  lambda <- D / Tt

  NV <- sum((lambda * a$tij - ped$ped_status)^2)
  W <- cumsum((lambda * a$T_j - a$D_j) / sqrt(NV))
  x <- cumsum(a$T_j) / Tt

  J <- length(x)
  cand <- which(x >= pi0 & x <= (1 - pi0) & seq_len(J) < J)
  if (!length(cand)) {
    return(list(stat = NA, p.value = 1, x = x, W = w))
  }

  w <- x[cand] * (1 - x[cand])
  LM <- W[cand]^2 / w
  stat <- max(LM)

  # Simulate reference borwnian bridges
  dx <- diff(c(0, x))
  Z <- matrix(rnorm(R * J), nrow = R, ncol = J) * rep(sqrt(dx), each = R)
  W <- t(apply(Z, 1, cumsum))
  B <- W - outer(W[, J], x)

  null <- apply(B[, cand, drop = FALSE]^2 / rep(w, each = R), 1, max)
  list(
    stat = stat,
    p.value = (1 + sum(null >= stat)) / (R + 1),
    x = x,
    W = w
  )
}

# Function to find split point
ll_node <- function(D, T_sum) {
  ifelse(D > 0, D * log(D / T_sum), 0) - D
}

best_split <- function(ped, minevents = 20) {
  a <- agg_ped(ped)
  Ds <- cumsum(a$D_j)
  Ts <- cumsum(a$T_j)
  D <- sum(a$D_j)
  Tt <- sum(a$T_j)
  J <- length(Ds)
  s <- seq_len(J - 1)
  ok <- Ds[s] >= minevents & (D - Ds[s]) >= minevents
  if (!any(ok)) {
    return(FALSE)
  }
  crit <- ll_node(Ds[s], Ts[s]) + ll_node(D - Ds[s], Tt - Ts[s])
  crit[!ok] <- -Inf
  a$tend[which.max(crit)]
}

pem_tree <- function(ped, alpha = 0.05, pi0 = 0.1, R = 1e4, minevents = 20) {
  D <- sum(ped$ped_status)
  leaf <- function(p) {
    data.frame(
      from = min(ped$tstart), to = max(ped$tend), events = D, lambda = D / sum(exp(ped$offset)), p.value = p
    )
  }
  if (D < 2 * minevents) {
    return(leaf(NA))
  }
  te <- suplm_time(ped, pi_0 = pi_0, R = R)
  if (is.na(te$stat) || te$p.value >= alpha) {
    return(leaf(te$p.value))
  }
  kappa <- best_split(ped, minevents = min_events)
  if (is.na(kappa)) {
    return(leaf(te$p.value))
  }
  rbind(
    pem_tree(ped[ped$tend <= kappa, ], alpha, pi0, R, minevents),
    pem_tree(ped[ped$tend > kappa, ], alpha, pi0, R, minevents)
  )
}

## Examples
# set.seed(1)
# beta0 <- -2
# delta <- 0.5
# tau <- 3
# sim_data <- sim_surv(
#   formula = ~ beta0 + delta * (t > tau),
#   data = data.frame(id = seq_len(1000)),
#   breaks = tau,
#   admin_time = 10,
#   lambda_c = 0.08
# )
# mean(sim_data$status == 0)
#
# ped1 <- as_ped(sim_data, Surv(time, status) ~ 1)
# res1_1 <- pem_tree(ped1)
# log(res1_1$lambda)
# source("R/pemtree_splines.R")
# res1_2 <- pemtree(
#   Surv(time, status) ~ offset(offset) | tend,
#   data = sim_data,
#   cut = NULL,
#   alpha = 0.05,
#   min_events = 20
# )
# res1_2$tree
# coef(res1_2$tree)
#
# set.seed(2)
# delta <- c(0.5, 0.7)
# tau <- c(2, 4)
# sim_data <- sim_surv(
#   formula = ~ beta0 + delta[1] * (t > tau[1]) + delta[2] * (t > tau[2]),
#   data = data.frame(id = seq_len(1000)),
#   breaks = tau,
#   admin_time = 10,
#   lambda_c = 0.05
# )
# mean(sim_data$status == 0)
#
# ped2 <- as_ped(sim_data, Surv(time, status) ~ 1)
# res2 <- pem_tree(ped2)
# res2$to
# log(res2$lambda)
