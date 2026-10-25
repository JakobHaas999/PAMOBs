# Shared model specification for the level-shift pilot and main studies.
make_level_shift_config <- function(phase, ncpus = 2L) {
  phase <- match.arg(phase, c("pilot", "main"))
  checkmate::assert_count(ncpus, positive = TRUE)

  list(
    run_id = paste0("level_shift_scenario_", phase, "_v1"),
    phase = phase,
    ncpus = ncpus,
    dgp = "generate_level_shift_data",
    fit_model = "fit_pemtree",
    evaluate = "evaluate_time_splits",
    packages = character(0),
    source = c(
      "R/sim_surv.R",
      "R/pemtree_splines.R",
      "simulation_batchtools/dgps/level_shift_hazard.R",
      "simulation_batchtools/methods/pemtree.R",
      "simulation_batchtools/metrics/tree_splits.R"
    ),
    prepare_pars = prepare_level_shift_pars
  )
}

prepare_level_shift_pars <- function(pars) {
  model_terms <- c(
    sprintf("splines::bs(tend, df = %d)", pars$df),
    if (pars$covariate_setting == "included") c("x1", "x2"),
    "offset(offset)"
  )
  # Only data columns and namespace-qualified functions are needed, plus
  # stats::offset. Do not serialize the whole study environment with each job.
  model_formula <- stats::as.formula(
    paste(
      "survival::Surv(time, status) ~",
      paste(model_terms, collapse = " + "),
      "| tend"
    ),
    env = asNamespace("stats")
  )

  # This grid is independent of the simulated change point.
  cuts <- sort(unique(c(
    seq(0, pars$admin_time, by = pars$cut_step),
    pars$admin_time
  )))

  pars$model_args <- list(
    formula = model_formula,
    cut = cuts,
    min_events = pars$min_events,
    alpha = pars$alpha,
    maxdepth = pars$maxdepth
  )
  pars
}
