config <- list(
  run_id = "pemtree_quadratic",
  phase = "main",
  ncpus = 2L
)
config$covariate_spec <- list(
  x1 = list(distfun = stats::rnorm),
  x2 = list(distfun = stats::rbinom, size = 1, prob = 0.5)
)
config$covariate_effects <- c(x1 = 0.3, x2 = 0.5)

config$model_args <- function(pars) {
  covariate_terms <- if (pars$covariate_setting == "included") {
    names(pars$covariate_spec)
  } else {
    character(0)
  }
  model_terms <- c("poly(tend, degree = 2, raw = TRUE)", covariate_terms, "offset(offset)")
  model_formula <- stats::as.formula(
    paste(
      "survival::Surv(time, status) ~",
      paste(model_terms, collapse = " + "),
      "| tend"
    )
  )

  list(formula = model_formula, min_events = pars$min_events, alpha = pars$alpha, maxdepth = pars$maxdepth)
}
config$dgp <- "quadratic_baseline"
config$fit_model <- "fit_pemtree"
config$evaluate <- "evaluate_time_splits"
config$source <- c(
  "R/sim_pexponential.R",
  "R/sim_censoring.R",
  "R/simulation_pam.R",
  "R/pemtree_splines.R",
  "simulation_batchtools/dgps/quadratic_baseline.R",
  "simulation_batchtools/methods/pemtree.R",
  "simulation_batchtools/metrics/tree_splits.R"
)
config$packages <- character(0)
config$prepare_pars <- function(pars) {
  pars$covariate_spec <- config$covariate_spec
  pars$covariate_effects <- config$covariate_effects
  pars$model_args <- config$model_args(pars)

  pars
}


design <- expand.grid(
  n = 1000L,

  beta0 = -1.5,
  beta1 = 1.2,
  beta2 = -0.2,

  delta = 0,
  tau = 0,

  covariate_setting = c("none", "included"),

  censoring_rate = 0.3,
  admin_time = 6,

  min_events = 25L,
  alpha = 0.05,
  maxdepth = 2L,

  rep = seq_len(1000L),

  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design$seed <- 100000L + seq_len(nrow(design))