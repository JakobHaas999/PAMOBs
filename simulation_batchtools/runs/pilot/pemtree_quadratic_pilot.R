# Quadratic-baseline pilot run: design, data generation and execution.
# Rscript starts the study; source() only defines the entry function.
run_pemtree_quadratic_pilot <- function(action = c("run", "prepare")) {
  action <- match.arg(action)

  config <- list(
    run_id = "pemtree_quadratic_pilot",
    phase = "pilot",
    ncpus = 2L,
    dgp = "quadratic_baseline",
    fit_model = "fit_pemtree",
    evaluate = "evaluate_time_splits",
    packages = character(0),
    source = c(
      "R/sim_pexponential.R",
      "R/sim_censoring.R",
      "R/simulation_pam.R",
      "R/pemtree_splines.R",
      "simulation_batchtools/R/fit_pemtree.R",
      "simulation_batchtools/R/tree_splits.R"
    )
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
    model_terms <- c(
      "poly(tend, degree = 2, raw = TRUE)",
      covariate_terms,
      "offset(offset)"
    )
    model_formula <- stats::as.formula(
      paste(
        "survival::Surv(time, status) ~",
        paste(model_terms, collapse = " + "),
        "| tend"
      ),
      env = asNamespace("stats")
    )
    list(
      formula = model_formula,
      min_events = pars$min_events,
      alpha = pars$alpha,
      maxdepth = pars$maxdepth
    )
  }
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
    rep = seq_len(20L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  design$seed <- 1000L + seq_len(nrow(design))

  quadratic_baseline <- function(pars) {
    log_hazard <- function(t) {
      value <- pars$beta0 + pars$beta1 * t + pars$beta2 * t^2
      if (pars$tau != 0) {
        value <- value + pars$delta * (t >= pars$tau)
      }
      value
    }
    if (pars$covariate_setting == "none") {
      sim_formula <- ~ log_hazard(t)
      pars$covariate_spec <- NULL
    } else {
      covariate_effects <- pars$covariate_effects[names(pars$covariate_spec)]
      covariate_terms <- sprintf(
        "(%s) * %s", covariate_effects, names(covariate_effects)
      )
      sim_formula <- stats::reformulate(
        termlabels = c("log_hazard(t)", covariate_terms)
      )
    }
    sim_pam(
      n = pars$n,
      formula = sim_formula,
      covariate_spec = pars$covariate_spec,
      censoring_rate = pars$censoring_rate,
      admin_time = pars$admin_time,
      cut = seq(0, pars$admin_time, by = 0.01)
    )
  }

  if (action == "prepare") {
    return(list(design = design, config = config))
  }

  source("simulation_batchtools/R/create_registry.R", local = TRUE)
  source("simulation_batchtools/R/run_registry.R", local = TRUE)
  source("simulation_batchtools/R/execute_simulation_run.R", local = TRUE)
  execute_simulation_run(design, config, envir = environment())
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 1L) stop("Supply at most one action: run or prepare.")
  invisible(run_pemtree_quadratic_pilot(if (length(args)) args[[1L]] else "run"))
}
