make_null_design <- function(n_rep) {
  

  scenarios <- data.frame(
    scenario_id = 1:3,
    scenario = c("cubic", "weibull", "lognormal"),
    df = c(6L, 6L, 5L),
  
    beta0 = c(-2.8,   NA_real_, NA_real_),
    beta1 = c( 0.25,  NA_real_, NA_real_),
    beta2 = c(-0.025, NA_real_, NA_real_),
    beta3 = c( 0.001, NA_real_, NA_real_),

    # Weibull-Hazard
    shape = c(NA_real_, 1.5,  NA_real_),
    rate  = c(NA_real_, 0.05, NA_real_),

    # Lognormal-Hazard
    meanlog = c(NA_real_, NA_real_, log(3)),
    sdlog   = c(NA_real_, NA_real_, 0.7),

    stringsAsFactors = FALSE
  )

  settings <- expand.grid(
    scenario_id = scenarios$scenario_id,
    n = c(500, 1000),
    censoring_rate = c(0.3, 0.5, 0.7),
    covariate_setting = c("none", "included"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  settings$setting_id <- seq_len(nrow(settings))

  settings <- merge(
  settings, scenarios, by = "scenario_id", sort = FALSE
  )

  design <- settings[rep(seq_len(nrow(settings)), each = n_rep), , drop = FALSE]
  design$rep <- rep(seq_len(n_rep), times = nrow(settings))
  rownames(design) <- NULL

  design$seed <- seq_len(nrow(design))

  # Add common parameters in all designs
  design$admin_time <- 10
  design$cut_step <- 0.1

  design$min_events <- 20L
  design$alpha <- 0.05
  design$maxdepth <- 2

  design$gamma1 <- 0.3
  design$gamm2 <- 0.5
  design$delta <- 0
  design$tau <- NA_real_

  return(design)

}