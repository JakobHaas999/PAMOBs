run_time_changing_effect_pilot <- function() {
  ncpus <- 2
  n_rep <- 5
  run_id <- "time_changing_effect_pilot_v1"
  dfs <- 6

  n_calibration <- 10000
  calibration_seed <- 20000000L
  seed_offset <- 10000000L

  source(
    "simulation_batchtools/R/study_helpers.R",
    local = TRUE
  )
  ## Baseline scenarios -------------------------------------------
  scenarios <- data.frame(
    scenario_id = 1:2,
    scenario = c("cubic", "weibull"),
    # Cubic log-hazard
    beta0 = c(-2.8, NA_real_),
    beta1 = c(0.25, NA_real_),
    beta2 = c(-0.025, NA_real_),
    beta3 = c(0.001, NA_real_),
    # Weibull hazard
    shape = c(NA_real_, 1.5),
    rate = c(NA_real_, 0.05)
  )

  ## Simulation design ---------------------------------------------
  settings <- expand.grid(
    scenario_id = scenarios$scenario_id,
    delta = seq(0, 0.9, by = 0.3),
    censoring_rate = c(0.4, 0.6),
    n = c(250, 500, 1000),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  settings <- cbind(
    settings,
    scenarios[
      match(settings$scenario_id, scenarios$scenario_id),
      setdiff(names(scenarios), "scenario_id"),
      drop = FALSE
    ]
  )
  rownames(settings) <- NULL

  settings$covariate_setting <- "included"
  settings$gamma1 <- 0.4
  settings$gamma2 <- 0.6

  settings$tau <- 3
  settings$admin_time <- 10

  settings$cut_step <- 0.1
  settings$min_events <- 20
  settings$alpha <- 0.05
  settings$maxdepth <- 2

  design <- cross_df_repetitions(settings, dfs, n_rep, seed_offset)

  ## DGP -------------------------------------------------------------
  # x1 changes effect with tau
  # x2 has a constant effect
  #
  # log h(t | x1, x2) = log h0(t) + gamma1 * x1 + gamma2 * x2 + delta * x1 * I(t >= tau)

  make_input <- function(pars, n = pars$n) {
    data <- data.frame(
      id = seq_len(n),
      x1 = rnorm(n),
      x2 = runif(n, min = -1, max = 1)
    )

    log_h0 <- switch(pars$scenario,
      cubic = function(t) {
        pars$beta0 + pars$beta1 * t + pars$beta2 * t^2 + pars$beta3 * t^3
      },
      weibull = function(t) {
        log(pars$rate) + log(pars$shape) + (pars$shape - 1) * log(t)
      },
      stop("Unknown scenario: ", pars$scenario)
    )

    formula <- ~ log_h0(t) + pars$gamma1 * x1 + pars$gamma2 * x2 + pars$delta * x1 * (t >= pars$tau)
    list(data = data, formula = formula)
  }

  generate_data <- function(pars) {
    input <- make_input(pars)

    sim_surv(
      formula = input$formula,
      data = input$data,
      breaks = if (pars$delta == 0) numeric(0) else pars$tau,
      admin_time = pars$admin_time,
      lambda_c = pars$lambda_c
    )
  }
  ## Registry configuration ---------------------------------------
  # prepare_spline_pars builds:
  #
  # Surv(time, status) ~
  #   bs(tend, df = 6) + x1 + x2 + offset(offset) | tend
  config <- list(
    run_id = run_id, phase = "pilot", ncpus = ncpus,
    dgp = "generate_data", fit_model = "fit_pemtree",
    evaluate = "evaluate_time_splits", packages = character(0),
    source = c(
      "R/sim_surv.R",
      "R/pemtree_splines.R",
      "simulation_batchtools/R/fit_pemtree.R",
      "simulation_batchtools/R/tree_splits.R"
    ),
    prepare_pars = prepare_spline_pars
  )

  ## Calibrate censoring --------------------------------------
  source("R/sim_surv.R", local = TRUE)

  source(
    "simulation_batchtools/R/censoring_calibration.R",
    local = TRUE
  )

  design <- calibrate_design(
    design = design,
    make_input = make_input,
    n_calibration = n_calibration,
    seed = calibration_seed
  )

  source(
    "simulation_batchtools/R/create_registry.R",
    local = TRUE
  )
  source(
    "simulation_batchtools/R/run_registry.R",
    local = TRUE
  )
  source(
    "simulation_batchtools/R/execute_simulation_run.R",
    local = TRUE
  )
  execute_simulation_run(
    design = design,
    config = config,
    envir = environment()
  )
}
if (sys.nframe() == 0L) {
  run_time_changing_effect_pilot()
}
