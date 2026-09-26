# Fit a PEM tree for temporal instability
#
# Fits a piecewise exponential model tree using time as the
# partitioning variable.
fit_temporal_pemtree <- function(data, cut = seq(0, 8, by = 0.1), ...) {
  pemtree(
    formula = Surv(time, status) ~ 1 + offset(offset) | tend,
    data = data,
    cut = cut,
    ...
  )
}
