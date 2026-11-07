Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("simulation_batchtools/R/complex_trees_helpers.R")

# ----- 1. Design -----
run_id <- "complex_tree_v1"
n_rep <- 1000
admin_time <- 10
n_calibration <- 50000
tau <- 3
beta <- -1

design <- CJ(
  n = 1000,
  scenario = c("covar_time_chang", "covar_struc_chang"),
  dz = c(0.5, 1),
  dt = c(1, 1.5),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(beta = beta, admin_time = admin_time, tau = tau)]

# ----- 2. Calibration -----
set.seed(32049)
calib <- unique(design[, .(scenario, dz, dt, censoring_rate)])[
  , calibrate_one(n_calibration, scenario, beta, dz, dt, tau, admin_time, censoring_rate),
  by = c("scenario", "dz", "dt")
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("scenario", "dz", "dt", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry -----
reg_dir <- "simulation_batchtools/registries/main"
dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, run_id),
  seed = 20303,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "Formula", "checkmate", "sandwich", "data.table", "simsurv"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R", "simulation_batchtools/R/complex_trees_helpers.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("complex_tree", fun = function(data, job, n, scenario, beta, dz, dt,
                                          tau, admin_time, lambda_c, ...) {
  generate_data(n, scenario, beta, dz, dt, tau, admin_time, lambda_c)
})
addAlgorithm("model_fit", fun = function(data, job, instance, cut_step, alpha,
                                         maxdepth, min_events, ...) {
  fit_and_summarise_complex_tree(instance, cut_step, alpha, maxdepth, min_events)
})
addExperiments(
  prob.designs = list(complex_tree = design),
  algo.designs = list(model_fit = data.table(
    cut_step = 0.1, alpha = 0.05, maxdepth = 4, min_events = 20
  )),
  repls = n_rep
)

# ----- 4. Test -----
pars <- unwrap(getJobPars())
res1 <- testJob(pars[dz == 1 & dt == 1 & scenario == "covar_time_chang"][1, job.id], external = TRUE)
res2 <- testJob(pars[dz == 1 & dt == 1.5 & scenario == "covar_struc_chang"][1, job.id], external = TRUE)
res1[, !"splits"]
res2[, !"splits"]


# ----- 5. Run ------
ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 25)]
submitJobs(ids)
waitForJobs()
