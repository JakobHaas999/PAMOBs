fit_tree_comparison <- function(dat, pars) {
  start <- proc.time()[["elapsed"]]
  args <- c(list(data = dat), pars$model_args)
  fit <- switch(
    pars$method,
    pemtree = do.call(pemtree, c(args, list(min_events = pars$min_events))),
    PErMtree = {
      set.seed(pars$seed + 1000000L)
      do.call(PErMtree, c(args, list(n_perm = pars$n_perm)))
    },
    stop("Unknown method: ", pars$method)
  )
  list(fit = fit, elapsed_sec = unname(proc.time()[["elapsed"]]- start))
}

evaluate_tree_comparison <- function(model, dat, pars) {
  fit <- model$fit
  if (pars$method == "pemtree") {
    n_nodes <- length(partykit::nodeids(fit$tree, terminal = TRUE))
    root <- partykit::node_party(fit$tree)
    p <- partykit::info_node(root)$p.value
  } else {
    n_nodes <- fit$terminal_nodes
    p <- fit$p_adjusted[["z"]]
  }
  if (length(p) != 1 || !is.finite(p)) p <- NA_real_
  data.frame(
    split_detected = n_nodes > 1,
    n_terminal_nodes = n_nodes,
    root_p_value = unname(p),
    test_available = is.finite(p),
    elapsed_sec = model$elapsed_sec,
    observed_censoring = mean(dat$status == 0)
  )
}