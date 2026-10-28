source("simulation_batchtools/R/study_helpers.R")

run_tree_comparison_pilot <- function(action = c("run", "prepare", "calibrate")) {
  action <- match.arg(action)
  n_rep <- 20
  n_calibration <- 10000

  settings <- expand.grid(
    scenario = c("cubic", "weibull"),
    delta = c(0, 0.4, 0.8),
    censoring_rate = c(0.4, 0.6),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  settings$n <- 1000
  settings$beta <- 0.4
  settings$prob_z <- 0.5
  settings$admin_time <- 10
  settings$cut_step <- 0.1
  settings$alpha <- 0.05
  settings$maxdepth <- 2
  settings$min_events <- 20
  settings$n_perm <- 499

  design <- cross_df_repetitions(settings, dfs = 6, n_rep = n_rep, seed_offset = 500000)
  design <- design[rep(seq_len(nrow(design)), each = 2), , drop = FALSE]
  design$method <- rep(c("pemtree", "PErMtree"), length.out = nrow(design))
  design$lambda_c <- NA_real_
  rownames(design) <- NULL

  make_input <- function(pars, n = pars$n) {
    dat <- data.frame(
      id = seq_len(n),
      x = rnorm(n),
      z = rbinom(n, size = 1, prob = pars$prob_z)
    )
    log_h0 <- switch(
      pars$scenario,
      cubic = function(t) -2.8 + 0.25 * t - 0.025 * t^2 + 0.001 * t^3,
      weibull = function(t) log(0.05) + log(1.5) + 0.5 * log(t)
    )
    formula <- ~ log_h0(t) + pars$beta * x + pars$delta * x * z
    list(data = dat, formula = formula)
  }

  generate_data <- function(pars) {
    input <- make_input(pars)
    dat <- sim_surv(
      formula = input$formula, data = input$data,
      admin_time = pars$admin_time, lambda_c = pars$lambda_c
    )
    dat$z <- factor(dat$z, levels = c(0, 1))
    dat
  }

  prepare_pars <- function(pars) {
    fit_formula <- as.formula(
      sprintf("survival::Surv(time, status) ~splines::bs(tend, df = %s) + x + offset(offset) | z", pars$df),
      env = asNamespace("stats")
    )
    pars$model_args <- list(
      formula = fit_formula,
      cut = seq(0, pars$admin_time, by = pars$cut_step),
      alpha = pars$alpha, maxdepth = pars$maxdepth
    )
    pars
  }

  config <- list(
    run_id = "tree_comparison_pilot_v1", phase = "pilot", ncpus = 2L,
    dgp = "generate_data", fit_model = "fit_tree_comparison",
    evaluate = "evaluate_tree_comparison", packages = character(),
    source =  c('R/sim_surv.R', 'R/pemtree_splines.R', 'R/permtree.R',
               'simulation_batchtools/R/tree_comparison.R'),
    prepare_pars = prepare_pars
  )

  if (action == "prepare") {
    message(nrow(design), " jobs;", length(unique(design$seed)), " simulated data sets")
    return(
      list(
        design = design, config = config, make_input, generate_data = generate_data
      )
    )
  }

  source('R/sim_surv.R', local = TRUE)
  for (id in unique(design$dgp_id)) {
    rows <- which(design$dgp_id == id)
    pars <- as.list(design[rows[1L], , drop = FALSE])
    message('Calibrating population ', id, '/', length(unique(design$dgp_id)))
    set.seed(920000L + id)
    input <- make_input(pars, n = n_calibration)
    events <- sim_event_times(
      formula = input$formula, data = input$data,
      maxt = pars$admin_time,
      interval = c(1e-8, 1.01 * pars$admin_time)
    )
    design$lambda_c[rows] <- calibrate_censoring(
      event_times = events, censoring_rate = pars$censoring_rate,
      admin_time = pars$admin_time
    )
  }

   if (action == 'calibrate') {
    return(invisible(list(design = design, config = config)))
  }
  source('simulation_batchtools/R/create_registry.R', local = TRUE)
  source('simulation_batchtools/R/run_registry.R', local = TRUE)
  source('simulation_batchtools/R/execute_simulation_run.R', local = TRUE)
  execute_simulation_run(design, config, envir = environment())
}

if (sys.nframe() == 0L) {
  action <- commandArgs(trailingOnly = TRUE)
  if (length(action) > 1L) stop('Supply at most one action: prepare, calibrate or run.')
  invisible(run_tree_comparison_pilot(if (length(action)) action[[1L]] else 'run'))
}