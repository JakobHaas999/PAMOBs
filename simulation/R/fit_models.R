# Fit a PEM tree for temporal instability
#
# Fits a piecewise exponential model tree using time as the
# partitioning variable.
fit_temporal_pemtree <- function(data, cut, min_events = 10, ...) {
  pemtree(
    formula = Surv(time, status) ~ 1 + offset(offset) | tend,
    data = data,
    cut = cut,
    min_events = min_events,
    ...
  )
}
