### Change for more complicated trees
evaluate_time_splits <- function(model, dat, pars) {
  tree <- model$tree
  n_terminal_nodes <- length(partykit::nodeids(tree, terminal = TRUE))
  split_detected <- n_terminal_nodes > 1
  root <- partykit::node_party(tree)
  root_p_value <- partykit::info_node(root)$p.value

  if (length(root_p_value) != 1 || !is.finite(root_p_value)) {
    root_p_value <- NA_real_
  }

  if (split_detected) {
    tau_hat <- partykit::breaks_split(partykit::split_node(root))[1]
  } else {
    tau_hat <- NA_real_
  }

  obs_cens <- mean(dat$status == 0)

  data.frame(
    split_detected = split_detected,
    n_terminal_nodes = n_terminal_nodes,
    root_p_value = root_p_value,
    tau_hat = tau_hat,
    abs_tau_error = if (pars$delta == 0) NA_real_ else abs(tau_hat - pars$tau),
    observed_censoring = obs_cens
  )
}