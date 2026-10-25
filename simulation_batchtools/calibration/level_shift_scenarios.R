# Run once before the pilot and main studies. Both use this same calibration.
local({
  source("simulation_batchtools/designs/level_shift_scenarios_base.R", local = TRUE)
  source("R/sim_surv.R", local = TRUE)
  source("simulation_batchtools/dgps/level_shift_hazard.R", local = TRUE)
  source("simulation_batchtools/R/level_shift_calibration.R", local = TRUE)

  design <- make_level_shift_design(n_rep = 1L)
  calibration <- calibrate_level_shift(design, n_pilot = 50000L, seed = 1900000L)

  output_file <- paste0(
    "simulation_batchtools/calibration/calibrations_results/",
    "level_shift_scenarios.rds"
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(calibration, output_file)
  message("Saved calibration for ", nrow(calibration), " settings to ", output_file)
})
