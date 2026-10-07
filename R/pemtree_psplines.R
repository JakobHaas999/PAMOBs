pemtree_psplines <- function(
  formula,
  data,
  cut,
  min_events = 20,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)

  f <- Formula::Formula(formula)
  if (length(f)[[2]] != 2) {
    stop("A two part formula with partitioning variables is required")
  }

  node_formula <- formula(f, lhs = 1, rhs = 1)
}
