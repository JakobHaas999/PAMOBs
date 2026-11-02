Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("simulation_batchtools/R/pemtree_covariate_helpers.R")

run_id <- "pemtree_covariate_main_v1"
n_rep <- 1000
beta <- 0.6
admin_time <- 10
n_calibration <- 50000

# ---- 1. Design ----------------------------------------------------------
design <- CJ(
  scenario = c("cubic", "weibull"),
  split_type = c("binary", "numeric"),
  n = c(250, 500, 1000),
  delta = seq(0, 1, by = 0.25),
  censoring_rate = seq(0.3, 0.7, by = 0.2)
)
design[, `:=`(beta = beta, admin_time = admin_time)]

# ---- 2. Calibration -----------------------------------------------------
set.seed(111)
calib <- unique(design[, .(scenario, split_type, delta)])[
  , calibrate_one(scenario, split_type, delta, beta, admin_time, unique(sort(design$censoring_rate)), n_calibration),
  by = .(scenario, split_type, delta)
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("scenario", "split_type", "delta", "censoring_rate"))
stopifnot(nrow(design) == n_before, !anyNA(design$lambda_c))

# ---- 3. Registry --------------------------------------------------------
dir.create("simulation_batchtools/registries/main", recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 2026,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "simsurv",
    "Formula", "checkmate", "sandwich"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R",
    "simulation_batchtools/R/pemtree_covariate_helpers.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("covariate", fun = function(data, job, n, scenario, split_type, delta,
                                       beta, admin_time, lambda_c, ...) {
  generate_data(n, scenario, split_type, delta, beta, admin_time, lambda_c)
})
addAlgorithm("pemtree", fun = function(data, job, instance, bs_df, cut_step, alpha,
                                       maxdepth, min_events, ...) {
  fit_and_summarise(instance, bs_df, cut_step, alpha, maxdepth, min_events)
})

addExperiments(
  prob.designs = list(covariate = design),
  algo.designs = list(pemtree = data.table(
    bs_df = 6, cut_step = 0.1, alpha = 0.05, maxdepth = 2, min_events = 20
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[split_type == "numeric" & n == 1000 & delta == 1][1, job.id], external = TRUE)

ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 25)]
submitJobs(ids)
waitForJobs()
