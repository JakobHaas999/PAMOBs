# null study: run-specific settings, design and data generation.
# Rscript <this-file> [prepare|calibrate|run]; source() only defines the function.
run_null_scenario_main <- function(action = c("run", "prepare", "calibrate")) {
  action <- match.arg(action)

  # Edit these values before starting a new registry.
  dfs <- c(4L, 6L, 8L)
  n_rep <- 1000L
  ncpus <- 2L
  run_id <- "null_scenario_main_df_v2"
  n_calibration <- 50000L
  calibration_seed <- 900000L
  calibration_file <- "simulation_batchtools/results/calibration/null_scenarios.rds"

  source("simulation_batchtools/R/study_helpers.R", local = TRUE)

  scenarios <- data.frame(
    scenario_id = 1:3,
    scenario = c("cubic", "weibull", "lognormal"),
    beta0 = c(-2.8, NA_real_, NA_real_),
    beta1 = c(0.25, NA_real_, NA_real_),
    beta2 = c(-0.025, NA_real_, NA_real_),
    beta3 = c(0.001, NA_real_, NA_real_),
    shape = c(NA_real_, 1.5, NA_real_),
    rate = c(NA_real_, 0.05, NA_real_),
    meanlog = c(NA_real_, NA_real_, log(3)),
    sdlog = c(NA_real_, NA_real_, 0.7)
  )
  settings <- expand.grid(
    scenario_id = scenarios$scenario_id,
    n = c(500L, 1000L),
    censoring_rate = c(0.3, 0.5, 0.7),
    covariate_setting = c("none", "included"),
    delta = 0,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  settings <- cbind(
    settings,
    scenarios[match(settings$scenario_id, scenarios$scenario_id),
              setdiff(names(scenarios), "scenario_id"), drop = FALSE]
  )
  settings$tau <- NA_real_
  settings$admin_time <- 10
  settings$gamma1 <- 0.3
  settings$gamm2 <- 0.5
  settings$cut_step <- 0.1
  settings$min_events <- 20L
  settings$alpha <- 0.05
  settings$maxdepth <- 2L

  design <- cross_df_repetitions(
    settings, dfs = dfs, n_rep = n_rep, seed_offset = 100000L
  )
  config <- list(
    run_id = run_id, phase = "main", ncpus = ncpus,
    dgp = "generate_data", fit_model = "fit_pemtree",
    evaluate = "evaluate_time_splits", packages = character(0),
    source = c(
      "R/sim_surv.R", "R/pemtree_splines.R",
      "simulation_batchtools/R/fit_pemtree.R",
      "simulation_batchtools/R/tree_splits.R"
    ),
    prepare_pars = prepare_spline_pars
  )

  make_input <- function(pars, n = pars$n) {
    checkmate::assert_count(n, positive = TRUE)
    checkmate::assert_number(pars$delta, finite = TRUE)
    if (pars$delta != 0) {
      checkmate::assert_number(pars$tau, finite = TRUE, lower = 0)
      if (pars$tau <= 0 || pars$tau >= pars$admin_time) {
        stop("tau must lie strictly between zero and admin_time.")
      }
    }
    dat <- data.frame(id = seq_len(n))
    log_h0 <- switch(
      pars$scenario,
      cubic = function(t) pars$beta0 + pars$beta1 * t + pars$beta2 * t^2 + pars$beta3 * t^3,
      weibull = function(t) log(pars$rate) + log(pars$shape) + (pars$shape - 1) * log(t),
      lognormal = function(t) {
        stats::dlnorm(t, pars$meanlog, pars$sdlog, log = TRUE) -
          stats::plnorm(t, pars$meanlog, pars$sdlog, lower.tail = FALSE, log.p = TRUE)
      },
      stop("Unknown scenario: ", pars$scenario)
    )
    log_hazard <- function(t) {
      log_h0(t) + if (pars$delta == 0) 0 else pars$delta * (t >= pars$tau)
    }
    if (pars$covariate_setting == "included") {
      dat$x1 <- stats::rnorm(n)
      dat$x2 <- stats::rbinom(n, size = 1, prob = 0.5)
      formula <- ~ log_hazard(t) + pars$gamma1 * x1 + pars$gamm2 * x2
    } else if (pars$covariate_setting == "none") {
      formula <- ~ log_hazard(t)
    } else {
      stop("Unknown covariate setting: ", pars$covariate_setting)
    }
    list(data = dat, formula = formula)
  }

  generate_data <- function(pars) {
    checkmate::assert_number(pars$lambda_c, finite = TRUE, lower = 0)
    input <- make_input(pars)
    sim_surv(
      formula = input$formula, data = input$data,
      breaks = if (pars$delta == 0) numeric(0) else pars$tau,
      admin_time = pars$admin_time, lambda_c = pars$lambda_c
    )
  }

  if (action == "prepare") {
    message(nrow(design), " planned fits; ", length(unique(design$seed)),
            " simulated data sets; df = ", paste(dfs, collapse = ", "))
    return(list(design = design, config = config,
                make_input = make_input, generate_data = generate_data))
  }

  source("R/sim_surv.R", local = TRUE)
  source("simulation_batchtools/R/censoring_calibration.R", local = TRUE)
  design <- load_or_calibrate_censoring(
    design, make_input, calibration_file,
    n_pilot = n_calibration, seed = calibration_seed
  )
  if (action == "calibrate") return(invisible(list(design = design, config = config)))

  source("simulation_batchtools/R/create_registry.R", local = TRUE)
  source("simulation_batchtools/R/run_registry.R", local = TRUE)
  source("simulation_batchtools/R/execute_simulation_run.R", local = TRUE)
  execute_simulation_run(design, config, envir = environment())
}

if (sys.nframe() == 0L) {
  action <- commandArgs(trailingOnly = TRUE)
  if (length(action) > 1L) stop("Supply at most one action: prepare, calibrate or run.")
  invisible(run_null_scenario_main(if (length(action)) action else "run"))
}
