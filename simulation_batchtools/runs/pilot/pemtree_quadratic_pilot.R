run_env <- new.env(parent = globalenv())

design_file <- "simulation_batchtools/designs/pilot/pemtree_quadratic_pilot.R"

source_files <- c(
  "simulation_batchtools/R/create_registry.R",
  "simulation_batchtools/R/run_registry.R",
  design_file
)

for (source_file in source_files) sys.source(source_file, envir = run_env)

registry_dir <- file.path(
  "simulation_batchtools", "registries", run_env$config$phase, run_env$config$run_id
)

if (!dir.exists(registry_dir)) run_env$create_registry(design_file)

reg <- run_env$run_registry(registry_dir = registry_dir, ncpus = run_env$config$ncpus)