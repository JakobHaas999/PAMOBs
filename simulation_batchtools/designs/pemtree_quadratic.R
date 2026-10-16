config <- list(
  run_id = "pemtree_null_quadratic",
  ncpus = 2L
)
config$covariate_spec <- list(
  x1 = list(distfun = stats::rnorm),
  x2 = list(stats::rbinom, size = 1, prob = 0.5)
)
config$covariate_effects <- c(x1 = 0.3, x2 = 0.5)

design <- expand.grid(
  n = 1000L,

  beta0 = -2,
  beta1 = 1.1,
  beta2 = -0.8,

  delta = 0,
  tau = 0,

  covariate_setting = c("none", "included"),

  censoring_rate = 0.3,
  admin_time = 10,

  min_events = 25L,
  alpha = 0.05,
  maxdepth = 2L,

  rep = seq_len(20L),

  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design$seed <- 1000L + seq_len(nrow(design))
