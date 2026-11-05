### setup
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/split_stability.R")
source("simulation_batchtools/R/split_stability_helpers.R")

# ----- 1. Design -----
run_id <- "split_stability_examples_v1"
n <- 1000
tau <- 3
beta <- 0.6
admin_time <- 10
design <- rbind(
  CJ(setting = "covariate", delta = c(0.5, 1), censoring_rate = c(0.3, 0.5, 0.7)),
  CJ(setting = "time_changing_effect", delta = c(0.75, 1), censoring_rate = c(0.3, 0.5, 0.7)),
  CJ(setting = "structural_change", delta = c(1, 1.5), censoring_rate = c(0.3, 0.5, 0.7))
)
design[, `:=`(n = n, tau = tau, beta = beta, admin_time = admin_time)]

# ----- 2. Calibration -----
set.seed(611)
n_calibration <- 50000
calib <- unique(design[, .(setting, delta, censoring_rate)])[
  , calibrate_one(setting, n_calibration, delta, tau, beta, admin_time, censoring_rate),
  by = c("setting", "delta")
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("setting", "delta", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry -----
dir.create("simulation_batchtools/registries/main", recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 2020,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "Formula",
    "checkmate", "sandwich", "data.table"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R", "R/split_stability.R", "simulation_batchtools/R/split_stability_helpers.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 8)

addProblem("examples_split_stability", fun = function(data, job, setting, n, delta,
                                                      tau, beta, lambda_c, admin_time, ...) {
  generate_experiment(setting, n, delta, tau, beta, lambda_c, admin_time)
})
addAlgorithm("model_fit", fun = function(data, job, instance, cut_step, alpha, min_events,
                                         B, frac, ...) {
  fit_and_summarise(
    instance,
    setting = job$prob.pars$setting,
    cut_step = cut_step,
    admin_time = job$prob.pars$admin_time,
    alpha = alpha,
    min_events = min_events,
    B = B,
    frac = frac
  )
})
addExperiments(
  prob.designs = list(examples_split_stability = design),
  algo.designs = list(model_fit = data.table(
    cut_step = 0.1, alpha = 0.05, min_events = 20, B = 500, frac = 0.632
  )),
  repls = 1
)

pars <- unwrap(getJobPars())
testJob(pars[setting == "structural_change" & delta == 1.5][1, job.id], external = TRUE)
submitJobs()
waitForJobs()