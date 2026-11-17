### Comparison of Goodman et al (2011) and exposure-scaled PEM tree
### setup
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_exact.R")
source("R/goodman_pch.R")

# ----- 1. Design -----
run_id <- "goodman_es_comparison_v1"
n_rep <- 1000
beta0 <- -2
tau <- 3
admin_time <- 10
n_calibration <- 50000

design <- CJ(
  n = c(250, 500, 1000),
  delta = seq(0, 1, by = 0.25),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(beta0 = beta0, tau = tau, admin_time = admin_time)]

# ----- 2. Calibration ------
set.seed(3100)
calib <- unique(design[, .(delta, censoring_rate)])[,
  {
    df <- data.frame(id = seq_len(n_calibration))
    sim_formula <- ~ beta0 + delta * (t >= tau)
    sim_data <- sim_surv(
      formula = sim_formula,
      data = df,
      breaks = tau,
      admin_time = admin_time,
      lambda_c = 0
    )
    event_times <- data.frame(eventtime = sim_data$time)

    data.table(
      censoring_rate = censoring_rate,
      lambda_c = vapply(censoring_rate, function(rate) {
        calibrate_censoring(event_times, rate, admin_time)
      }, numeric(1))
    )
  },
  by = c("delta")
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("delta", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry -----
reg_dir <- "simulation_batchtools/registries/main"
dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, run_id),
  seed = 20203,
  packages = c(
    "survival", "pammtools", "simsurv",
    "checkmate", "data.table", "eventTrack"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_exact.R", "R/goodman_pch.R"
  )
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("piecewise_constant_hazard", fun = function(data, job, n, beta0, delta, tau,
                                                       admin_time, lambda_c, ...) {
  df <- data.frame(id = seq_len(n))
  sim_surv(
    formula = ~ beta0 + delta * (t >= tau),
    data = df,
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
})

addAlgorithm("models_fit", fun = function(data, job, instance, alpha, maxdepth, min_events,
                                          pi0, R, K_max, goodman_levels, ...) {
  data <- instance[, c("time", "status")]

  fits <- vector("list", length = 2)
  fits[[1]] <- tryCatch(
    {
      ped <- as_ped(data, Surv(time, status) ~ 1)
      pem_tree(
        ped,
        alpha = alpha, pi0 = pi0, R = R, minevents = min_events, maxdepth = maxdepth
      )
    },
    error = function(e) e
  )
  fits[[2]] <- tryCatch(
    {
      goodman_forward(
        time = data$time, status = data$status,
        K_max = K_max, alpha = alpha, levels = goodman_levels
      )
    },
    error = function(e) e
  )

  summarise_pem_tree <- function(fit) {
    if (inherits(fit, "error")) {
      return(data.table(
        model = "pemtree_es", split = NA, p_root = NA_real_,
        stat = NA_real_, split_point = NA_real_,
        lambda_1 = NA_real_, lambda_2 = NA_real_,
        error = conditionMessage(fit)
      ))
    }

    root <- fit[fit$depth == 1, ]
    split <- !root$terminal
    p <- root$p.value
    data.table(
      model = "pemtree_es", split = split,
      p_root = if (length(p) == 1) p else NA_real_,
      stat = NA_real_,
      split_point = if (split) root$split else NA_real_,
      lambda_1 = if (split) fit[fit$depth == 2, ]$lambda[1] else NA_real_,
      lambda_2 = if (split) fit[fit$depth == 2, ]$lambda[2] else NA_real_,
      error = NA_character_
    )
  }

  summarise_goodman <- function(fit) {
    if (inherits(fit, "error") || !isTRUE(fit$ok) || is.null(fit$tests)) {
      msg <- if (inherits(fit, "error")) conditionMessage(fit) else "goodman fit failed"
      return(data.table(
        model = "goodman", split = NA, p_root = NA_real_,
        stat = NA_real_, split_point = NA_real_,
        lambda_1 = NA_real_, lambda_2 = NA_real_,
        error = msg
      ))
    }
    split <- fit$k_hat >= 1
    data.table(
      model = "goodman", split = split, p_root = fit$tests$p[1],
      stat = fit$tests$Xw[1],
      split_point = if (split) fit$tau[1] else NA_real_,
      lambda_1 = if (split) fit$lambda[1] else NA_real_,
      lambda_2 = if (split) fit$lambda[2] else NA_real_,
      error = NA_character_
    )
  }

  rbind(summarise_pem_tree(fits[[1]]), summarise_goodman(fits[[2]]))[
    , observed_censoring := mean(data$status == 0)
  ][]
})

addExperiments(
  prob.designs = list(piecewise_constant_hazard = design),
  algo.designs = list(models_fit = data.table(
    alpha = 0.05, maxdepth = 2, min_events = 20,
    pi0 = 0.1, R = 1e4, K_max = 1, goodman_levels = "goodman"
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[n == 1000 & delta == 1 & censoring_rate == 0.3][1, job.id], external = TRUE)
ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 50)]
submitJobs(ids)
waitForJobs()
