design_file <- "simulation_batchtools/designs/pilot/null_scenarios_pilot.R"
source(design_file)
source("simulation_batchtools/R/create_registry.R")
source("simulation_batchtools/R/run_registry.R")

registry_dir <- file.path("simulation_batchtools", "registries", config$phase, config$run_id)

if (!dir.exists(registry_dir)) {
  create_registry(design_file)
}

reg <- run_registry(registry_dir = registry_dir, ncpus = config$ncpus)
