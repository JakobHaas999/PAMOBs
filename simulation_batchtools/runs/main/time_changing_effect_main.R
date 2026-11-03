Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("simulation_batchtools/R/time_changing_effect_helpers.R")


# ----- 1. Design ------
run_id <- "time_changing_effect_main_v2"
n_rep <- 1000
admin_time <- 10
n_calibration <- 50000
tau <- 3
beta <- 0.6

design <- CJ(
  n = c(250, 500, 1000),
  scenario = c("cubic", "weibull"),
  delta = seq(0, 1, by = 0.25),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(beta = beta, admin_time = admin_time, tau = tau)]

# ------ 2. Calibration -----
set.seed(11111)
calib <- unique(design[, .(scenario, delta, censoring_rate)])[
  , calibrate_one(scenario, delta, tau, beta, admin_time, censoring_rate, n_calibration),
  by = c("scenario", "delta")
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("scenario", "delta", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry -----
dir.create("simulation_batchtools/registries/main", recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 3000,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "simsurv", "Formula", "checkmate", "sandwich", "data.table"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R", "simulation_batchtools/R/time_changing_effect_helpers.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("time_changing_effect", fun = function(data, job, n, scenario, delta, tau, beta,
                                                  admin_time, lambda_c, ...) {
  generate_data(n, scenario, delta, tau, beta, admin_time, lambda_c)
})
addAlgorithm("fit_model", fun = function(data, job, instance, bs_df, cut_step,
                                         alpha, maxdepth, min_events, ...) {
  fit_and_summarise(instance, bs_df, cut_step, alpha, maxdepth, min_events)
})
addExperiments(
  prob.designs = list(time_changing_effect = design),
  algo.designs = list(fit_model = data.table(
    bs_df = 6, cut_step = 0.1, alpha = 0.05, maxdepth = 2, min_events = 20
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[n == 1000 & delta == 1][1, job.id], external = TRUE)
submitJobs()
waitForJobs()
