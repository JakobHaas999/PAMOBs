create_registry <- function(design_file) {
  study_env <- new.env(parent = globalenv())
  sys.source(design_file, envir = study_env)

  config <- study_env$config
  design <- study_env$design

  if (!isTRUE(config$phase %in% c("pilot", "main"))) {
    stop("config$phase must be 'pilot' or 'main'.")
  }

  source_files <- c("simulation_batchtools/R/run_one.R", config$source)
  for (source_file in source_files) sys.source(source_file, envir = study_env)
  
  get_fun <- function(name) get(name, envir = study_env, mode = "function", inherits = FALSE)
  runner <- get_fun("run_one")
  dgp <- get_fun(config$dgp)
  fit_model <- get_fun(config$fit_model)
  evaluate <- get_fun(config$evaluate)

  parameter_sets <- lapply(seq_len(nrow(design)), function(i) {
    pars <- as.list(design[i, , drop = FALSE])
    config$prepare_pars(pars)
  })

  registry_root <- file.path("simulation_batchtools", "registries", config$phase)
  registry_dir <- file.path(registry_root, config$run_id)
  if (file.exists(registry_dir)) {
    stop("Registry already exists: ", registry_dir, ".",
    "\nLoad existing registry or take a new run_id")
  }

  dir.create(registry_root, recursive = TRUE, showWarnings = FALSE)

  # Save registry
  reg <- batchtools::makeRegistry(
    file.dir = registry_dir,
    work.dir = getwd(),
    conf.file = NA,
    packages = config$packages,
    source = character(0),
    make.default = FALSE
  )

  ids <- batchtools::batchMap(
    fun = runner,
    pars = parameter_sets,
    more.args = list(dgp = dgp, fit_model = fit_model, evaluate = evaluate),
    reg = reg
  )

  job_design <- cbind(as.data.frame(ids), design)

  saveRDS(job_design, file.path(reg$file.dir, "job_design.rds"))
  message(nrow(ids), " Jobs registered in: ", reg$file.dir)
  invisible(reg)
}
