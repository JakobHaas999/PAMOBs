source("simulation_batchtools/designs/null_scenarios_base.R")
design <- make_null_design(n_rep = 1L)
source("R/sim_surv.R")
source("simulation_batchtools/dgps/null_hazard.R")

n_pilot <- 50000L

groups <- split(design, interaction(design$scenario_id, design$covariate_setting, drop = TRUE))

calibration <- lapply(seq_along(groups), function(i) {
  group <- groups[[i]]
  pars <- as.list(group[1, , drop = FALSE])
  
  calibration_seed <- 900000L + i
  set.seed(calibration_seed)

  input <- make_null_input(pars, n = n_pilot)
  events <- sim_event_times(
    formula = input$formula, data = input$data,
    maxt = pars$admin_time, interval = c(1e-8, 1.01 * pars$admin_time)
  )

  targets <- unique(group[c("setting_id", "censoring_rate")])
  targets$lambda_c <- vapply(targets$censoring_rate, function(target) {
    calibrate_censoring(event_times = events, censoring_rate = target, admin_time = pars$admin_time)
  }, numeric(1))

  targets$calibration_seed <- calibration_seed

  targets
})

calibration <- do.call(rbind, calibration)
rownames(calibration) <- NULL

saveRDS(calibration, "simulation_batchtools/calibration/calibrations_results/null_scenarios.rds")