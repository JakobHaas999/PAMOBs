construct_tree_formula <- function(covariate_spec = NULL, baseline_terms = character()) {
  model_terms <- c(
    baseline_terms, names(covariate_spec), "offset(offset)"
  )

  as.formula(
    paste("Surv(time, status) ~", paste(model_terms, collapse = " + "), "| tend"),
    env = parent.frame()
  )
}

extract_time_split <- function(fit) {
  root <- partykit::node_party(fit$tree)
  split <- partykit::split_node(root)

  if (is.null(split)) {
    return(NA_real_)
  }

  partykit::breaks_split(split)[[1]]
}
