### setup
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/split_stability.R")

# ----- 1. Design -----
run_id <- "split_stability_piecewise_constant_v1"
n <- 1000
tau <- 3
beta0 <- -2
admin_time <- 10
n_calibration <- 50000
design <- CJ(
  delta = c(0.5, 0.75, 1),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(beta0 = beta0, admin_time = admin_time, tau = tau, n = n)]

# ----- 2. Calibration
set.seed(30302)
calib <- unique(design[, .(delta, censoring_rate)])[,
  {
    data <- data.frame(id = seq_len(n_calibration))
    sim_formula <- ~ beta0 + delta * (t >= tau)
    d <- sim_surv(
      sim_formula, data,
      breaks = tau, admin_time = admin_time, lambda_c = 0
    )
    event_times <- data.frame(eventtime = d$time)
    data.table(
      censoring_rate = censoring_rate,
      lambda_c = vapply(censoring_rate, function(rate) {
        calibrate_censoring(event_times, rate, admin_time)
      }, numeric(1))
    )
  },
  by = "delta"
]
n_before <- nrow(design)
merge(design, calib, by = c("delta", "censoring_rate"))
stopifnot(nrow(design) == n_before)

# ----- 3. Registry -----
dir.create("simulation_batchtools/registries/main", recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path("simulation_batchtools/registries/main", run_id),
  seed = 23213,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "Formula",
    "checkmate", "sandwich", "data.table"
  ),
  source = c("R/sim_surv.R", "R/pemtree_splines.R", "R/split_stability.R")
)
reg$cluster.functions <- makeClusterFunctionsSocket(ncpus = 2)

addProblem("piecewise_constant_hazard", fun = function(data, job, n, delta, censoring_rate, beta0,
                                                       admin_time, tau, lambda_c, ...) {
  data <- data.frame(id = seq_len(n))
  sim_formula <- ~ beta0 + delta * (t >= tau)
  sim_surv(
    formula = sim_formula,
    data = data,
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
})
