config <- list(
  run_id = "pemtree_null_quadratic",
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

  rep = seq_len(20L),

  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design$seed <- 1000L + seq_len(nrow(design))
