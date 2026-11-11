wald_pch <- function(time, status, K = 3, alpha = 0.05,
                     scheme = c("goodman", "eventTrack"),
                     refit = TRUE) {
  scheme <- match.arg(scheme)
  a <- if (scheme == "goodman") 2 * alpha else alpha

  fit <- tryCatch(eventTrack::piecewiseExp_MLE(time, status, K = K), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(ok = FALSE, k_hat = NA_integer_))
  }

  tab <- eventTrack::piecewiseExp_test_changepoint(fit, alpha = a)
  rej <- tab$reject == 1
  k_hat <- sum(cumprod(rej))

  if (k_hat == 0) {
    tau <- numeric(0)
    lambda <- sum(status) / sum(time)
  } else if (refit && k_hat < K) {
    f2 <- tryCatch(eventTrack::piecewiseExp_MLE(time, status, k_hat), error = function(e) NULL)
    if (is.null(f2)) {
      return(list(ok = FALSE, k_hat = NA_integer_))
    }
    tau <- f2$tau
    lambda <- f2$lambda
  } else {
    tau <- fit$tau[seq_len(k_hat)]
    lambda <- fit$lambda[seq_len(k_hat + 1)]
  }

  list(
    ok = TRUE, k_hat = k_hat, hit_K = k_hat == K, tau = tau,
    lambda = lambda, Xw = tab$Xw, p = tab[["p-value"]]
  )
}

goodman_forward <- function(time, status, K_max, alpha = 0.05,
                            levels = c("goodman", "eventTrack")) {
  levels <- match.arg(levels)

  out <- list(
    ok = TRUE, fail_at = NA_integer_, k_hat = 0, hit_K = FALSE,
    tau = numeric(0), lambda = sum(status) / sum(time), tests = NULL
  )

  for (k in seq_len(K_max)) {
    fit <- tryCatch(eventTrack::piecewiseExp_MLE(time, status, K = k), error = function(e) NULL)
    if (is.null(fit)) {
      out$ok <- FALSE
      out$fail_at <- k
      break
    }

    lam <- fit$lambda
    se <- fit$lambda.SE
    Xw <- (head(lam, -1) - tail(lam, -1))^2 / (head(se, -1)^2 + tail(se, -1)^2)

    if (any(!is.finite(Xw))) {
      out$ok <- FALSE
      out$fail_at <- k
      break
    }

    j <- which.min(Xw)
    stat <- Xw[[j]]
    p <- pchisq(stat, df = 1, lower.tail = FALSE)
    level <- if (levels == "goodman") alpha / 2^(k - 1) else alpha / 2^k

    if (!is.finite(p)) {
      out$ok <- FALSE
      out$fail_at <- k
      break
    }

    out$tests <- rbind(
      out$tests,
      data.frame(k = k, pair = j, Xw = stat, p = p, level = level, reject = p < level)
    )

    if (p >= level) break

    out$k_hat <- k
    out$tau <- fit$tau
    out$lambda <- fit$lambda
  }

  out$hit_K <- out$k_hat == K_max
  out
}

## Example
# ## 1) One break point
# # Simulate data
# set.seed(1)
# n <- 1000
# beta0 <- -2.5
# delta <- 0.5
# tau <- 5
# admin_time <- 10
#
# source("R/sim_surv.R")
# sim_data <- sim_surv(
#   formula = ~ beta0 + delta * (t > tau),
#   data = data.frame(id = seq_len(n)),
#   breaks = tau,
#   lambda_c = 0.06,
#   admin_time = admin_time
# )
# # Check censoring
# mean(sim_data$status == 0)
# #> [1] 0.49
#
# # Goodman et al (2011) method
# res_goodman <- goodman_forward(
#   time = sim_data$time, status = sim_data$status,
#   K_max = 5, alpha = 0.05,
#   levels = "goodman"
# )
#
# # Pemtree
# source("R/pemtree_splines.R")
# res_pemtree <- pemtree(
#   formula = survival::Surv(time, status) ~ offset(offset) | tend,
#   data = sim_data,
#   cut = seq(0, ceiling(max(sim_data$time)), by = 0.1),
#   min_events = 20,
#   alpha = 0.05,
#   maxdepth = 5
# )
#
# # Compare results
# tree <- res_pemtree$tree
# k_hat_pem <- partykit::width(tree) - 1
# inner_nodes <- setdiff(partykit::nodeids(tree), partykit::nodeids(tree, terminal = TRUE))
# tau_pem <- sort(unlist(partykit::nodeapply(tree,
#   ids = inner_nodes,
#   FUN = function(node) partykit::breaks_split(partykit::split_node(node))
# )))
# c(goodman = res_goodman$k_hat, pemtree = k_hat_pem)
# rbind(goodman = res_goodman$tau, pemtree = tau_pem)
#
# ## 2) Two break points
# # Simulate data
# set.seed(2)
# n <- 1000
# beta0 <- -2.5
# delta <- c(0.5, 0.6)
# tau <- c(3, 6)
# admin_time <- 10
#
# sim_data <- sim_surv(
#   formula = ~ beta0 + delta[1] * (t > tau[1]) + delta[2] * (t > tau[2]),
#   data = data.frame(id = seq_len(n)),
#   breaks = tau,
#   lambda_c = 0.06,
#   admin_time = admin_time
# )
# # Check censoring
# mean(sim_data$status == 0)
# #> [1] 0.381
#
# # Goodman et al (2011) method
# res_goodman <- goodman_forward(
#   time = sim_data$time, status = sim_data$status,
#   K_max = 5, alpha = 0.05,
#   levels = "goodman"
# )
#
# # Pemtree
# res_pemtree <- pemtree(
#   formula = Surv(time, status) ~ offset(offset) | tend,
#   data = sim_data,
#   cut = seq(0, ceiling(max(sim_data$time)), by = 0.1),
#   min_events = 20,
#   alpha = 0.05,
#   maxdepth = 5
# )
#
# # Compare results
# tree <- res_pemtree$tree
# k_hat_pem <- partykit::width(tree) - 1
# inner_nodes <- setdiff(nodeids(tree), nodeids(tree, terminal = TRUE))
# tau_pem <- sort(unlist(partykit::nodeapply(
#   tree,
#   ids = inner_nodes,
#   FUN = function(node) partykit::breaks_split(partykit::split_node(node))
# )))
#
# c(goodman = res_goodman$k_hat, pemtree = k_hat_pem)
# max_tau_length <- max(length(res_goodman$tau), length(tau_pem))
# tau_g <- res_goodman$tau
# length(tau_g) <- length(tau_pem) <- max_tau_length
# res_tau <- rbind(goodman = tau_g, pemtree = tau_pem)
# colnames(res_tau) <- NULL
# res_tau
