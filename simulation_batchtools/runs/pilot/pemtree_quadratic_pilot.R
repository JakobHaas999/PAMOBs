# Run from the project root.
registry_env <- new.env(parent = globalenv())

sys.source(
  "simulation_batchtools/R/create_registry.R",
  envir = registry_env
)

reg <- registry_env$create_registry(
  design_file = "simulation_batchtools/designs/pilot/pemtree_quadratic_pilot.R"
)

batchtools::getStatus(reg = reg)
