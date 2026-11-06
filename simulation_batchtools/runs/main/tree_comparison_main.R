Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("simulation_batchtools/R/pemtree_covariate_helpers.R")
source("simulation_batchtools/R/tree_comparison_helpers.R")

# ----- 1. Design -----
run_id <- "tree_comparison_v1"
n_rep <- 1000
beta <- 0.6
admin_time <- 10
n_calibration <- 50000

design <- CJ(
  scenario = c("cubic", "weibull"),
  split_type = c("binary", "numeric"),
  n = c(250, 500, 1000),
  delta = seq(0, 1, by = 0.25),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(beta = beta, admin_time = admin_time)]

# ----- 2. Calibration -----
set.seed(1633)
calib <- unique(design[, .(scenario, split_type, delta, censoring_rate)])[
  , calibrate_one(scenario, split_type, delta, beta, admin_time, censoring_rate, n_calibration),
  by = c("scenario", "split_type", "delta")
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("scenario", "split_type", "delta", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry -----
dir.create("simulation_batchtools/registries/main", showWarnings = FALSE, recursive = TRUE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 1634,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "simsurv",
    "Formula", "checkmate", "sandwich", "data.table"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R",
    "simulation_batchtools/R/pemtree_covariate_helpers.R",
    "simulation_batchtools/R/tree_comparison_helpers.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("covariate", fun = function(data, job, n, scenario, split_type, delta,
                                       beta, admin_time, lambda_c, ...) {
  generate_data(n, scenario, split_type, delta, beta, admin_time, lambda_c)
})
addAlgorithm("model_fit", fun = function(data, job, instance, bs_df, cut_step, alpha,
                                         maxdepth, min_events, ...) {
  fit_and_summarise_tree_comp(instance, bs_df, cut_step, alpha, maxdepth, min_events)
})
addExperiments(
  prob.designs = list(covariate = design),
  algo.designs = list(model_fit = data.table(
    bs_df = 6, cut_step = 0.1, alpha = 0.05, maxdepth = 2, min_events = 20
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[n == 1000 & delta == 0.75 & split_type == "numeric"][1, job.id], external = TRUE)

ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 25)]
submitJobs(ids)
waitForJobs()
