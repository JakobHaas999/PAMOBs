### Final run: subsampling stability of the split point
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/pemtree_exact.R")
source("R/split_stability.R")

# ----- 1. Design -----
run_id <- "split_stability_final_v1"
n <- 1000
tau <- 3
beta0 <- -2
beta <- 0.6
admin_time <- 10
n_calibration <- 50000

design <- CJ(
  setting = c("covariate", "pch"),
  delta = c(0.5, 0.75),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(n = n, tau = tau, beta0 = beta0, beta = beta, admin_time = admin_time)]

# ----- 2. Calibration -----
set.seed(123432)
calib <- unique(design[, .(setting, delta, censoring_rate)])[,
  {
    if (setting == "covariate") {
      z <- round(runif(n_calibration, min = -1, max = 1), 2)
      df <- data.frame(
        id = seq_len(n_calibration),
        x = rnorm(n_calibration),
        z = z, g = as.integer(z > 0)
      )
      sim_formula <- ~ log(1.5) + log(0.05) + 0.5 * log(t) + beta * x + delta * g * x
      sim_data <- sim_surv(
        formula = sim_formula,
        data = df,
        admin_time = admin_time,
        lambda_c = 0
      )
    } else {
      df <- data.frame(id = seq_len(n_calibration))
      sim_formula <- ~ beta0 + delta * (t >= tau)
      sim_data <- sim_surv(
        formula = sim_formula,
        data = df,
        breaks = tau,
        admin_time = admin_time,
        lambda_c = 0
      )
    }

    event_times <- data.frame(eventtime = sim_data$time)
    data.table(
      censoring_rate = censoring_rate,
      lambda_c = vapply(censoring_rate, function(r) {
        calibrate_censoring(event_times, r, admin_time = admin_time)
      }, numeric(1))
    )
  },
  by = c("setting", "delta")
]
design <- merge(design, calib, by = c("setting", "delta", "censoring_rate"))

# ----- 3. Registry -----
reg_path <- "simulation_batchtools/registries/main"
dir.create(reg_path, recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_path, run_id),
  seed = 67676,
  packages = c(
    "survival", "partykit", "pammtools", "simsurv", "splines", "Formula",
    "checkmate", "sandwich", "data.table"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R", "R/split_stability.R", "R/pemtree_exact.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("split_stability", fun = function(data, job, setting, n, delta, tau, beta0,
                                             beta, lambda_c, admin_time, ...) {
  if (setting == "covariate") {
    z <- round(runif(n, min = -1, max = 1), digits = 2)
    df <- data.frame(
      id = seq_len(n),
      x = rnorm(n),
      z = z,
      g = as.integer(z > 0)
    )
    sim_formula <- ~ log(1.5) + log(0.05) + 0.5 * log(t) + beta * x + delta * x * g
    return(sim_surv(
      formula = sim_formula, data = df,
      admin_time = admin_time, lambda_c = lambda_c
    ))
  } else {
    df <- data.frame(id = seq_len(n))
    sim_formula <- ~ beta0 + delta * (t >= tau)
    return(sim_surv(
      formula = sim_formula, data = df, breaks = tau,
      admin_time = admin_time, lambda_c = lambda_c
    ))
  }
})

addAlgorithm("stability", fun = function(data, job, instance, alpha, B, min_events,
                                         frac, pi0, R, keep_all, ...) {
  pp <- job$prob.pars
  res <- if (pp$setting == "pch") {
    stability_pem_tree(
      instance[, c("time", "status")],
      cut = seq(0, 10, by = 0.1), B = B, frac = frac,
      alpha = alpha, pi0 = pi0, R = R, min_events = min_events
    )
  } else {
    fit_args <- list(
      formula = Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | z,
      data = instance[, c("time", "status", "x", "z")],
      cut = seq(0, 10, by = 0.1),
      alpha = alpha,
      min_events = min_events
    )
    res <- split_stability(pemtree, fit_args = fit_args, B = B, frac = frac)
    root <- node_party(res$full_fit$tree)
    full_split <- !is.terminal(root)
    data.table(
      b = seq_len(B),
      split = res$split,
      split_point = res$split_point,
      full_split = full_split,
      full_split_point = if (full_split) breaks_split(split_node(root)) else NA_real_
    )
  }
  res[, observed_censoring := mean(instance$status == 0)]
  if (keep_all) {
    return(res)
  }

  truth <- if (pp$setting == "covariate") 0 else pp$tau
  sp <- res[split == TRUE, split_point]
  q <- if (length(sp) > 0) quantile(sp, c(0.025, 0.975), names = FALSE) else c(NA_real_, NA_real_)
  data.table(
    split_rate = mean(res$split, na.rm = TRUE),
    n_split = length(sp),
    full_split = res$full_split[[1]],
    full_split_point = res$full_split_point[[1]],
    lower = q[[1]], upper = q[[2]], width = q[[2]] - q[[1]],
    covered = truth >= q[[1]] & truth <= q[[2]],
    observed_censoring = res$observed_censoring[[1]]
  )
})

algo <- function(B, keep_all) {
  data.table(
    alpha = 0.05, min_events = 20, B = B, frac = 0.632,
    pi0 = 0.1, R = 1e4, keep_all = keep_all
  )
}

# Part A: Illustration
addExperiments(
  prob.designs = list(split_stability = design),
  algo.designs = list(stability = algo(B = 500, keep_all = TRUE)),
  repls = 1
)
# Part B: Coverage
addExperiments(
  prob.designs = list(split_stability = design[censoring_rate == 0.5]),
  algo.designs = list(stability = algo(B = 200, keep_all = FALSE)),
  repls = 500
)

pars <- unwrap(getJobPars())
testJob(pars[setting == "pch" & keep_all == FALSE][1, job.id], external = TRUE)
testJob(pars[setting == "covariate"][1, job.id], external = TRUE)

ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 20)]
submitJobs(ids)
waitForJobs()
