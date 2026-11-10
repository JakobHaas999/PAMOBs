#### Small simulation to detect splitpoints in a constant hazard
### setup
library("batchtools")
library("data.table")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

# ----- 1. Design -----
run_id <- "pemtree_piecewise_constant_hazard_v1"
n_rep <- 1000
beta0 <- -2.5
admin_time <- 10
n_calibration <- 50000
tau <- 5

design <- CJ(
  n = 500,
  delta = c(0, 0.5, 1),
  censoring_rate = 0.5
)
design[, `:=`(beta0 = beta0, admin_time = admin_time, tau = tau)]

# ----- 2. Calibration -----
set.seed(4094)
calib <- unique(design[, .(delta, censoring_rate)])[,
  {
    data <- data.frame(id = seq_len(n_calibration))
    surv_data <- sim_surv(
      formula = ~ beta0 + delta * (t >= tau),
      data = data,
      breaks = tau,
      admin_time = admin_time,
      lambda_c = 0
    )
    event_times <- data.frame(eventtime = surv_data$time)
    lambda_c <- calibrate_censoring(event_times, censoring_rate, admin_time = admin_time)
    data.table(censoring_rate = censoring_rate, lambda_c = lambda_c)
  },
  by = "delta"
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("delta", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry -----
dir.create("simulation_batchtools/registries/main", recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 393920,
  packages = c(
    "survival", "partykit", "pammtools", "simsurv", "Formula", "checkmate",
    "sandwich", "data.table"
  ),
  source = c("R/sim_surv.R", "R/pemtree_splines.R")
)
reg$cluster.functions <- makeClusterFunctionsSocket(ncpus = 2)

addProblem("constant_hazard", fun = function(data, job, n, beta0, delta, tau,
                                             admin_time, lambda_c, ...) {
  data <- data.frame(id = seq_len(n))
  sim_surv(
    formula = ~ beta0 + delta * (t >= tau),
    data = data,
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
})
addAlgorithm("model_fit", fun = function(data, job, instance, cut_step, alpha,
                                         maxdepth, min_events, ...) {
  data <- instance[, c("time", "status")]
  fit <- tryCatch(
    pemtree(
      formula = Surv(time, status) ~ offset(offset) | tend,
      data = data,
      min_events = min_events,
      cut = seq(0, ceiling(max(data$time)), by = cut_step),
      alpha = alpha,
      maxdepth = maxdepth
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.table(
      split = NA, p_root = NA_real_, split_point = NA_real_,
      error = conditionMessage(fit), observed_censoring = mean(data$status == 0)
    ))
  }

  root <- node_party(fit$tree)
  split <- !is.terminal(root)
  p <- info_node(root)$p.value
  split_point <- if (split) breaks_split(split_node(root)) else NA_real_

  data.table(
    split = split, p_root = if (length(p) == 1) p else NA_real_,
    split_point = split_point, error = NA_character_,
    observed_censoring = mean(data$status == 0)
  )
})

addExperiments(
  prob.designs = list(constant_hazard = design),
  algo.designs = list(model_fit = data.table(
    cut_step = 0.1, alpha = 0.05, maxdepth = 2, min_events = 20
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
for (d in unique(design$delta)) {
  print(testJob(pars[delta == d][1, job.id], external = TRUE))
}

submitJobs()
waitForJobs()
