pamtree <- function(
  formula,
  data,
  cut = NULL,
  min_events,
  gam_method = "REML",
  sp,
  cluster_id = NULL,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)

  # Asserts
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_numeric(cut, lower = 0, any.missing = FALSE)
  checkmate::assert_count(min_events, positive = TRUE)
  checkmate::assert_string(cluster_id, null.ok = TRUE)

  # Construct PED formula
  f <- Formula::Formula(formula)
  if (length(f)[[2]] != 2) {
    stop("A two part formula is partitioning variables is required")
  }

  node_formula <- formula(f, lhs = 1, rhs = 1)
  node_formula <- update(node_formula, ped_status ~ .)


  data_vars <- intersect(all.vars(f), colnames(data))
  response <- f[[2]]
  ped_vars <- setdiff(data_vars, all.vars(response))
  ped_formula <- reformulate(
    termlabels = ped_vars,
    response = deparse(response)
  )

  ped_args <- list(formula = ped_formula, data = data)
  ped <- pammtools::as_ped(ped_formula, data = data, cut = cut)
}
