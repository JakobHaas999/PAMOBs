run_covariate_scenarios <- function() {
  n_rep <- 20
  n_calibration <- 10000

  settings <- expand.grid(
    scenario = c("cubic", "weibull"),
    split_variable = c("binary", "numeric"),
    delta = c(0, 0.3, 0.6, 0.9),
    n = c(500, 1000),
    censoring_rate = c(0.3, 0.5, 0.7),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  settings$beta <- 0.6
  settings$prob_binary <- 0.5
  settings$admin_time <- 10
  settings$cut_step <- 0.1
  settings$alpha <- 0.05
  settings$maxdepth <- 2
  settings$min_events <- 20

  


}