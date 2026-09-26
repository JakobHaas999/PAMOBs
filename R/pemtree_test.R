pemtree <- function(
  formula,
  data,
  cut,
  node_formula = ~1,
  partition = ~tend,
  ...
) {
  ped <- pammtools::as_ped(
    formula = formula,
    data = data,
    cut = cut
  )

  # Extract rhs expressions
  node_rhs <- attr(te)
}
