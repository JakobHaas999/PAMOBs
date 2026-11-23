Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
### Run: Multiple split points (2 to 3) in a piecewise constant hazard
### setup
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_exact.R")
source("R/goodman_pch.R")

# ----- 1. Design -----
run_id <- "multiple_split_points_v1"
n_rep <- 1000
beta0 <- -2
admin_time <- 10
n_calibration <- 50000

mcp_scenario <- function(scenario) {
  switch(scenario,
    stair = list(tau = c(1.5, 3.5), delta = c(0.5, 0.5)),
    up_and_down = list(tau = c(1.5, 3.5), delta = c(0.5, -0.5)),
    three_splits = list(tau = c(1, 2.5, 4), delta = c(0.6, -0.4, 0.5)),
    stop("Unknown scenario")
  )
}

mcp_formula <- function(beta0, tau, delta) {
  as.formula(
    paste("~", paste(c(beta0, sprintf("(%s) * (t > %s)", delta, tau)), collapse = " + "))
  )
}

design <- CJ(
  n = c(500, 1000),
  censoring_rate = 0.5,
  scenario = c("stair", "up_and_down", "three_splits")
)
design[, `:=`(beta0 = beta0, admin_time = admin_time)]

# ----- 2. Calibration -----
set.seed(23232)
calib <- unique(design[, .(scenario, censoring_rate)])[,
  {
    sc <- mcp_scenario(scenario)
    sim_data <- sim_surv(
      formula = mcp_formula(beta0, sc$tau, sc$delta),
      data.frame(id = seq_len(n_calibration)),
      breaks = sc$tau, admin_time = admin_time, lambda_c = 0
    )
    .(lambda_c = calibrate_censoring(
      data.frame(eventtime = sim_data$time), censoring_rate, admin_time
    ))
  },
  by = c("scenario", "censoring_rate")
]
design <- merge(design, calib, by = c("scenario", "censoring_rate"))

# ----- 3. Registry -----
reg_dir <- "simulation_batchtools/registries/main"
dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, run_id),
  seed = 20205,
  packages = c(
    "survival", "pammtools", "simsurv", "checkmate", "data.table", "eventTrack"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_exact.R", "R/goodman_pch.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("multiple_change_points", fun = function(data, job, scenario, n, beta0,
                                                    admin_time, lambda_c, ...) {
  sc <- switch(scenario,
    stair = list(tau = c(1.5, 3.5), delta = c(0.5, 0.5)),
    up_and_down = list(tau = c(1.5, 3.5), delta = c(0.5, -0.5)),
    three_splits = list(tau = c(1, 2.5, 4), delta = c(0.6, -0.4, 0.5)),
    stop("Unknown scenario")
  )
  f <- as.formula(
    paste("~", paste(c(beta0, sprintf("(%s) * (t > %s)", sc$delta, sc$tau)), collapse = " + "))
  )
  sim_surv(
    formula = f,
    data = data.frame(id = seq_len(n)),
    breaks = sc$tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
})

addAlgorithm("models_fit", fun = function(data, job, instance, alpha, maxdepth,
                                          min_events, pi0, R, K_max, ...) {
  data <- instance[, c("time", "status")]
  out <- function(model, splits = numeric(0), error = NA_character_) {
    data.table(
      model = model,
      n_splits = if (is.na(error)) length(splits) else NA_integer_,
      split_points = list(sort(splits)),
      error = error
    )
  }

  tree <- tryCatch(
    {
      fit <- pem_tree(
        ped = as_ped(Surv(time, status) ~ 1, data = data),
        alpha = alpha, pi0 = pi0, R = R,
        minevents = min_events, maxdepth = maxdepth
      )
      out("pemtree_es", fit$split[!fit$terminal])
    },
    error = function(e) out("pemtree_es", error = conditionMessage(e))
  )

  goodman <- tryCatch(
    {
      fit <- goodman_forward(
        time = data$time, status = data$status,
        K_max = K_max, alpha = alpha,
        levels = "goodman"
      )
      if (!isTRUE(fit$ok) && fit$k_hat == 0) stop("goodman fit failed")
      out("goodman", splits = fit$tau)
    },
    error = function(e) out("goodman", error = conditionMessage(e))
  )

  rbind(tree, goodman)[, observed_censoring := mean(data$status == 0)][]
})

addExperiments(
  prob.designs = list(multiple_change_points = design),
  algo.designs = list(models_fit = data.table(
    alpha = 0.05, maxdepth = Inf, min_events = 20, pi0 = 0.1, R = 1e4, K_max = 4
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[scenario == "three_splits" & n == 1000][1, job.id], external = TRUE)
ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 50)]
submitJobs(ids)
waitForJobs()
