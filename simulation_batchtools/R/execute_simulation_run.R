# An existing registry may only be resumed with exactly the same job design.
execute_simulation_run <- function(design, config, envir = parent.frame()) {
  registry_dir <- file.path("simulation_batchtools", "registries", config$phase, config$run_id)
  if (dir.exists(registry_dir)) {
    design_file <- file.path(registry_dir, "job_design.rds")
    if (!file.exists(design_file)) stop("Existing registry has no saved job design.")
    saved <- as.data.frame(readRDS(design_file))
    saved$job.id <- NULL
    if (!isTRUE(all.equal(saved, design, check.attributes = FALSE, tolerance = 0))) {
      stop("Existing registry has a different design. Choose a new config$run_id.")
    }
  } else {
    create_registry(design, config, envir = envir)
  }
  run_registry(registry_dir = registry_dir, ncpus = config$ncpus)
}
