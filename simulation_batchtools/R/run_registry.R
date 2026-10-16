run_registry <- function(registry_dir, ncpus) {
  reg <- batchtools::loadRegistry(
    file.dir = registry_dir, writeable = TRUE, conf.file = NA, make.default = FALSE
  )
  ids <- batchtools::findNotSubmitted(reg = reg)
  if (!nrow(ids)) {
    message("No new jobs to start exists.")
    batchtools::getStatus(reg = reg)
    return(invisible(reg))
  }

  reg$cluster.functions <- batchtools::makeClusterFunctionsSocket(ncpus = ncpus)
  
  batchtools::submitJobs(ids = ids, reg = reg)
  completed <- batchtools::waitForJobs(ids = ids, reg = reg)
  batchtools::getStatus(reg = reg)

  if (!isTRUE(completed)) {
    warning("Not all submiited jobs completed successfully.", call. = FALSE)
  }

  invisible(reg)
}
