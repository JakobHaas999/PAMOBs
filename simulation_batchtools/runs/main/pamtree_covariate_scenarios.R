### setup
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/pamtree.R")
source("simulation_batchtools/R/pamtree_covariate_helpers.R")

# Fit pemtree, pamtree, weibull AFT tree to covariate scenarios

# ----- 1. Design -----
run_id <- "pamtree_covariate_scenarios_v1"
n_rep <- 1000
beta <- 0.6
admin_time <- 10
n_calibration <- 50000

design <- CJ(
  scenario = c("cubic", "weibull"),
  split_type = c("binary", "numeric"),
  n = c(250, 500),
  delta = c(0, 0.5, 1),
  censoring_rate = seq(0.3, 0.7, by = 0.2)
)
design[, `:=`(beta = beta, admin_time = admin_time)]

# ----- 2. Calibration -----
set.seed(20202)
calib <- unique(design[, .(scenario, split_type, delta, censoring_rate)])[
  , calibrate_one(scenario, split_type, delta, beta, admin_time, censoring_rate, n_calibration),
  by = c("scenario", "split_type", "delta")
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("scenario", "split_type", "delta", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry -----
dir.create("simulation_batchtools/registries/main", recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 4930,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "simsurv",
    "Formula", "checkmate", "sandwich", "data.table", "mgcv"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R", "R/pamtree.R",
    "simulation_batchtools/R/pemtree_covariate_helpers.R",
    "simulation_batchtools/R/tree_comparison_helpers.R",
    "simulation_batchtools/R/pamtree_covariate_helpers.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("covariate", fun = function(data, job, n, scenario, split_type, delta,
                                       beta, admin_time, lambda_c, ...) {
  generate_data(n, scenario, split_type, delta, beta, admin_time, lambda_c)
})

addAlgorithm("model_fit", fun = function(data, job, instance, bs_df, gam_method, sp, cut_step,
                                         alpha, maxdepth, min_events, ...) {
  fit_and_summarise_pam(instance, bs_df, gam_method, sp, cut_step, alpha, maxdepth, min_events)
})

addExperiments(
  prob.designs = list(covariate = design),
  algo.designs = list(model_fit = data.table(
    bs_df = 6, gam_method = "REML", sp = NA_real_, cut_step = 0.1, alpha = 0.05, maxdepth = 2, min_events = 20
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[n == 500 & delta == 1 & split_type == "numeric"][1, job.id], external = TRUE)

ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 20)]
submitJobs(ids)
waitForJobs()
