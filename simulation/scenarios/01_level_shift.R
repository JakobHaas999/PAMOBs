# Scenario 01: Level shift in baseline hazard
#
# The log-baseline hazard is piecewise constant with a
# single structural change at tau:
#
# log(lambda(t)) = beta0 + delta * I(t >= tau)
#
# delta = 0 corresponds to the null scenario without parameter instability

scenario_level_shift <- function(
  n,
  delta,
  tau,
  beta0,
  censoring_rate = NULL,
  covariate_spec = NULL,
  covariate_effects = NULL,
  admin_time = 10,
  sim_interval = 0.05
) {
  if (is.null(admin_time)) {
    stop("admin_rate must be specified")
  }
  if (tau > admin_time) {
    stop("tau must not exceed admin_time")
  }

  if (is.null(covariate_spec)) {
    if (!is.null(covariate_effects)) {
      stop("covariate_effects must be NULL when covariate_spec is NULL")
    }
    covariate_terms <- character(0)
  } else {
    covariate_names <- names(covariate_spec)
    if (is.null(covariate_names) ||
      is.null(names(covariate_effects)) ||
      !setequal(covariate_names, names(covariate_effects))) {
      stop("covariate_effects must be a named numeric vector matching covariate_spec")
    }

    covariate_effects <- covariate_effects[covariate_names]
    covariate_terms <- sprintf("(%s) * %s", covariate_effects, covariate_names)
  }

  log_hazard <- function(t) {
    beta0 + delta * (t >= tau)
  }

  sim_formula <- reformulate(
    termlabels = c("log_hazard(t)", covariate_terms),
    env = environment()
  )

  sim_cut <- sort(unique(c(
    seq(0, admin_time, by = sim_interval),
    tau,
    admin_time
  )))

  sim_pam(
    n = n,
    formula = sim_formula,
    covariate_spec = covariate_spec,
    censoring_rate = censoring_rate,
    admin_time = admin_time,
    cut = sim_cut
  )
}
