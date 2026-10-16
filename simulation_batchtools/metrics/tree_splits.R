### Change for more complicated trees
evaluate_time_splits <- function(model, dat, pars) {
  tree <- model$tree
  n_terminal_nodes <- length(partykit::nodeids(tree, terminal = TRUE))
  split_detected <- n_terminal_nodes > 1
  if (split_detected) {
    root <- partykit::node_party(tree)
    tau_hat <- partykit::breaks_split(partykit::split_node(root))[1]
  } else {
    tau_hat <- NA_real_
  }

  obs_cens <- mean(dat$status == 0)

  data.frame(
    split_detected = split_detected,
    n_terminal_nodes = n_terminal_nodes,
    tau_hat = tau_hat,
    abs_tau_error = if (pars$delta) NA_real_ else abs(tau_hat - pars$tau),
    observed_censoring = obs_cens
  )
}