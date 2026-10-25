source("simulation_batchtools/designs/level_shift_scenarios_base.R", local = TRUE)
source("simulation_batchtools/R/level_shift_calibration.R", local = TRUE)
source("simulation_batchtools/R/level_shift_config.R", local = TRUE)

design <- make_level_shift_design(n_rep = 1000L)
calibration_file <- paste0(
  "simulation_batchtools/calibration/calibrations_results/",
  "level_shift_scenarios.rds"
)
design <- attach_level_shift_calibration(design, calibration_file)

config <- make_level_shift_config(phase = "main", ncpus = 2L)
