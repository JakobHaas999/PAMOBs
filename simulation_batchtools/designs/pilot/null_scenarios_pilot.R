source("simulation_batchtools/designs/null_scenarios_base.R")
design <- make_null_design(n_rep = 5L)


calibration_file <- "simulation_batchtools/calibration/calibrations_results/null_scenarios.rds"
if (file.exists("simulation_batchtools/calibration/calibrations_results/null_scenarios.rds")) {
  calibration <- readRDS(calibration_file)
} else {
  source("simulation_batchtools/calibration/null_scenarios.R")
  calibration <- readRDS(calibration_file)
}

index <- match(design$setting_id, calibration$setting_id)
design$lambda_c <- calibration$lambda_c[index]
stopifnot(
  "lambda_c" %in% names(design),
  all(is.finite(design$lambda_c)),
  all(design$lambda_c >= 0)
)

config <- list(
  "run_id" = "null_scenario_pilot_v1",
  phase = "pilot",
  ncpus = 2L,
  dgp = "generate_null_data",
  fit_model = "fit_pemtree",
  evaluate = "evaluate_time_splits",

  packages = character(0),

  source = c(
    "R/sim_surv.R",
    "R/pemtree_splines.R",
    "simulation_batchtools/dgps/null_hazard.R",
    "simulation_batchtools/methods/pemtree.R",
    "simulation_batchtools/metrics/tree_splits.R"

  )
)

config$prepare_pars <- function(pars) {
  model_terms <- c(
    sprintf("splines::bs(tend, df = %d)", pars$df),
    if (pars$covariate_setting == "included") c("x1", "x2"),
    "offset(offset)"
  )
  model_formula <- stats::as.formula(
    paste(
      "survival::Surv(time, status) ~", paste(model_terms, collapse = " + "), "| tend"
    )
  )
  cuts <- sort(unique(c(
    seq(0, pars$admin_time, by = pars$cut_step),
    pars$admin_time
  )))

  pars$model_args <- list(
    formula = model_formula, cut = cuts,
    min_events = pars$min_events, alpha = pars$alpha, maxdepth = pars$maxdepth
  )

  pars
}

