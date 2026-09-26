# Scenario 01: Level shift in baseline hazard
#
# The log-baseline hazard is piecewise constant with a
# single structural change at tau:
#
# log(lambda(t)) = beta0 + delta * I(t > tau)
#
# delta = 0 corresponds to the null scenario without parameter instability

scenario_level_shift <- function(
  n,
  delta,
  tau,
  beta0,
  censoring_rate = NULL
) {
  log_hazard <- function(t) {
    beta0 + delta * ifelse(t > tau, delta, 0)
  }

  sim_pam(
    n = n,
    formula = ~ log_hazard(t),
    covariate_spec = NULL,
    censoring_rate = censoring_rate,
    cut = seq(0, 8, by = 0.1)
  )
}
