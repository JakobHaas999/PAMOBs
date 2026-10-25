# Run from the project root. Small real simulations use a temporary registry;
# no production calibration, registry or result file is created or changed.
local({
  source("simulation_batchtools/designs/level_shift_scenarios_base.R", local = TRUE)
  source("simulation_batchtools/R/level_shift_config.R", local = TRUE)
  source("simulation_batchtools/R/level_shift_calibration.R", local = TRUE)
  source("simulation_batchtools/dgps/null_hazard.R", local = TRUE)
  source("simulation_batchtools/dgps/level_shift_hazard.R", local = TRUE)
  source("R/sim_surv.R", local = TRUE)

  settings <- make_level_shift_design(n_rep = 1L)
  pilot <- make_level_shift_design(n_rep = 100L)
  main <- make_level_shift_design(n_rep = 1000L)
  stopifnot(
    nrow(settings) == 144L,
    nrow(pilot) == 14400L,
    nrow(main) == 144000L,
    !anyDuplicated(main$seed),
    all(pilot$seed %in% main$seed),
    all(table(main$setting_id) == 1000L)
  )

  # Check the actual log-hazard formula, including both covariate effects.
  population_rows <- !duplicated(settings[c("scenario", "covariate_setting", "delta")])
  for (i in which(population_rows)) {
    pars <- as.list(settings[i, , drop = FALSE])
    input <- make_level_shift_input(pars, n = 1L)
    reference <- make_null_input(pars, n = 1L)
    values <- list(t = c(4, 5, 6), x1 = 0.25, x2 = 1)
    eta <- eval(input$formula[[2L]], values, environment(input$formula))
    eta_null <- eval(reference$formula[[2L]], values, environment(reference$formula))
    stopifnot(isTRUE(all.equal(
      eta - eta_null, pars$delta * (values$t >= pars$tau), tolerance = 1e-12
    )))
  }

  # Only six populations and 100 calibration subjects per population here.
  # These small calibration samples must never be used for production runs.
  small <- settings[
    settings$n == 500 & settings$censoring_rate == 0.7 &
      settings$delta == log(3),
    , drop = FALSE
  ]
  small$n <- 200L
  calibration <- calibrate_level_shift(small, n_pilot = 100L, seed = 1800000L)
  calibration_file <- tempfile("level-shift-smoke-calibration-", fileext = ".rds")
  saveRDS(calibration, calibration_file)
  calibrated <- attach_level_shift_calibration(small, calibration_file)

  # A setting ID alone must not allow calibration for a different true shift.
  changed <- small
  changed$delta <- log(2)
  stale <- try(attach_level_shift_calibration(changed, calibration_file), silent = TRUE)
  stopifnot(inherits(stale, "try-error"))

  config <- make_level_shift_config("pilot")
  job_env <- new.env(parent = globalenv())
  for (file in c("simulation_batchtools/R/run_one.R", config$source)) {
    sys.source(file, envir = job_env)
  }
  pars <- lapply(seq_len(nrow(calibrated)), function(i) {
    config$prepare_pars(as.list(calibrated[i, , drop = FALSE]))
  })
  reg <- batchtools::makeRegistry(
    file.dir = tempfile("level-shift-smoke-registry-"),
    work.dir = getwd(),
    conf.file = NA,
    make.default = FALSE
  )
  ids <- batchtools::batchMap(
    fun = job_env$run_one,
    pars = pars,
    more.args = list(
      dgp = job_env$generate_level_shift_data,
      fit_model = job_env$fit_pemtree,
      evaluate = job_env$evaluate_time_splits
    ),
    reg = reg
  )

  results <- lapply(ids$job.id, function(id) {
    result <- batchtools::testJob(id, external = TRUE, reg = reg)
    if (!isTRUE(result$success)) stop("Smoke job failed: ", result$error)
    stopifnot(
      nrow(result) == 1L,
      !is.na(result$split_detected),
      result$n_terminal_nodes %in% c(1L, 2L),
      is.finite(result$observed_censoring),
      result$observed_censoring >= 0,
      result$observed_censoring <= 1
    )
    result
  })
  message(
    "Level-shift checks passed: 144 settings, 24 hazard specifications, ",
    "calibration validation and six independent R-process fits."
  )
})
