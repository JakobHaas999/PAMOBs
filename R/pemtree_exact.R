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

suplm_time <- function(ped, pi0 = 0.1, R = 1e4) {
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
    return(list(stat = NA, p.value = 1, x = x, W = W))
  }

  w <- x[cand] * (1 - x[cand])
  LM <- W[cand]^2 / w
  stat <- max(LM)

  # Simulate reference Brownian bridges (J x R)
  dx <- diff(c(0, x))
  Z <- matrix(rnorm(J * R), nrow = J, ncol = R) * sqrt(dx)
  W_sim <- apply(Z, 2, cumsum)
  B_sim <- W_sim - outer(x, W_sim[J, ])
  null <- apply(B_sim[cand, , drop = FALSE]^2 / w, 2, max)

  list(
    stat = stat,
    p.value = (1 + sum(null >= stat)) / (R + 1),
    x = x,
    W = W
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
    return(NA)
  }
  crit <- ll_node(Ds[s], Ts[s]) + ll_node(D - Ds[s], Tt - Ts[s])
  crit[!ok] <- -Inf
  a$tend[which.max(crit)]
}

pem_tree <- function(ped, alpha = 0.05, pi0 = 0.1, R = 1e4, minevents = 20, maxdepth = Inf, depth = 1) {
  D <- sum(ped$ped_status)
  node <- function(p, split = NA_real_) {
    data.frame(
      depth = depth,
      terminal = is.na(split),
      from = min(ped$tstart),
      to = max(ped$tend),
      events = D,
      lambda = D / sum(exp(ped$offset)),
      p.value = p,
      split = split
    )
  }
  if (D < 2 * minevents || depth >= maxdepth) {
    return(node(NA))
  }
  te <- suplm_time(ped, pi0 = pi0, R = R)
  if (is.na(te$stat) || te$p.value >= alpha) {
    return(node(te$p.value))
  }
  kappa <- best_split(ped, minevents = minevents)
  if (is.na(kappa)) {
    return(node(te$p.value))
  }
  rbind(
    node(te$p.value, kappa),
    pem_tree(ped[ped$tend <= kappa, ], alpha, pi0, R, minevents, maxdepth, depth + 1),
    pem_tree(ped[ped$tend > kappa, ], alpha, pi0, R, minevents, maxdepth, depth + 1)
  )
}
