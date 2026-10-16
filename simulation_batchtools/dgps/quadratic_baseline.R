quadratic_baseline <- function(pars) {

  log_hazard <- function(t) {
    fun <- pars$beta0 + pars$beta1 * t + pars$beta2 * t^2
    if (pars$tau != 0) {
      fun <- fun + pars$delta * (t >= pars$tau)
    }
    fun
  }

  if (pars$covariate_setting == "none") {
    sim_formula <- ~ log_hazard(t)
    pars$covariate_spec <- NULL
  } else {
    covar_effects <- pars$covariate_effects[names(pars$covariate_spec)]
    covar_terms <- sprintf("(%s) * %s", covar_effects, names(covar_effects))
    sim_formula <- reformulate(
      termlabels = c("log_hazard(t)", covar_terms)
    )
  }

  sim_pam(
    n = pars$n,
    formula = sim_formula,
    covariate_spec = covariate_spec,
    censoring_rate = pars$censoring_rate,
    admin_time = pars$admin_time,
    cut = seq(0, pars$admin_time, by = 0.01)
  )
}