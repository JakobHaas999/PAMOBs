Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/pamtree.R")
source("simulation_batchtools/R/pamtree_null_helpers.R")

# ----- 1. Design -----
run_id <- "pamtree_null_event_grid_v1"
n_rep <- 200
admin_time <- 10
n_calibration <- 50000
n_obs <- 1000
censoring_rate <- 0.5

design <- data.table(
  censoring_rate = censoring_rate, admin_time = admin_time, n = n_obs, scenario = "weibull_hazard"
)

# ----- 2. Calibration -----
set.seed(22323)
calib <- calibrate_one(n_calibration, admin_time, censoring_rate)
design$lambda_c <- calib$lambda_c

# ------ 3. Registry -----
dir.create("simulation_batchtools/registries/main", recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 23323,
  packages = c("survival", "partykit", "pammtools", "simsurv", "Formula", "checkmate", "sandwich", "mgcv", "data.table"),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R", "R/pamtree.R", "simulation_batchtools/R/pamtree_null_helpers.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("weibull_hazard", fun = function(data, job, n, admin_time, lambda_c, ...) {
  generate_data(n, admin_time, lambda_c)
})
addAlgorithm("model_fit", fun = function(data, job, instance, cut_step, alpha, maxdepth, min_events, ...) {
  fit_and_summarise(instance, cut_step, alpha, maxdepth, min_events)
})
addExperiments(
  prob.designs = list(weibull_hazard = design),
  algo.designs = list(model_fit = data.table(
    cut_step = NA_real_, alpha = 0.05, maxdepth = 2, min_events = 20
  )),
  repls = n_rep
)

testJob(1, external = TRUE)
submitJobs()
waitForJobs()
