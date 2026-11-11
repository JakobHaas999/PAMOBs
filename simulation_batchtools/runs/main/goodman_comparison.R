### Pemtree comparison with Goodman et al (2011)
### setup
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/goodman_pch.R")

# ----- 1. Design -----
run_id <- "goodman_comparison_v1"
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
set.seed(30404)
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
  by = "delta"
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("delta", "censoring_rate"))
stopifnot(n_before == nrow(design))

# ----- 3. Registry ------
reg_path <- "simulation_batchtools/registries/main"
dir.create(reg_path, recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_path, run_id),
  seed = 40404,
  packages = c(
    "survival", "partykit", "pammtools", "simsurv", "Formula", "checkmate",
    "sandwich", "data.table", "eventTrack"
  ),
  source = c("R/sim_surv.R", "R/pemtree_splines.R", "R/goodman_pch.R")
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

addAlgorithm("models_fit", fun = function(data, job, instance, cut_step, alpha, maxdepth,
                                          min_events, K_max, goodman_levels, ...) {
  data <- instance[, c("time", "status")]

  fits <- vector("list", length = 2)
  fits[[1]] <- tryCatch(
    pemtree(
      Surv(time, status) ~ offset(offset) | tend,
      data = data,
      cut = seq(0, ceiling(max(data$time)), by = cut_step),
      min_events = min_events,
      alpha = alpha,
      maxdepth = maxdepth
    ),
    error = function(e) e
  )

  fits[[2]] <- tryCatch(
    goodman_forward(
      time = data$time, status = data$status,
      K_max = K_max, alpha = alpha, levels = goodman_levels
    ),
    error = function(e) e
  )

  summarise_pemtree <- function(fit) {
    if (inherits(fit, "error")) {
      return(data.table(
        model = "pemtree", split = NA, p_root = NA_real_, stat = NA_real_,
        split_point = NA_real_, error = conditionMessage(fit)
      ))
    }
    root <- node_party(fit$tree)
    split <- !is.terminal(root)
    p <- info_node(root)$p.value
    data.table(
      model = "pemtree", split = split,
      p_root = if (length(p) == 1) unname(p) else NA_real_,
      stat = NA_real_,
      split_point = if (split) breaks_split(split_node(root)) else NA_real_,
      error = NA_character_
    )
  }

  summarise_goodman <- function(fit) {
    if (inherits(fit, "error") || !isTRUE(fit$ok) || is.null(fit$tests)) {
      msg <- if (inherits(fit, "error")) conditionMessage(fit) else "goodman fit failed"
      return(data.table(
        model = "goodman", split = NA, p_root = NA_real_, stat = NA_real_,
        split_point = NA_real_, error = msg
      ))
    }
    split <- fit$k_hat >= 1
    data.table(
      model = "goodman", split = split, p_root = fit$tests$p[1], stat = fit$tests$Xw[1],
      split_point = if (split) fit$tau[1] else NA_real_, error = NA_character_
    )
  }
  rbind(summarise_pemtree(fits[[1]]), summarise_goodman(fits[[2]]))[
    , observed_censoring := mean(data$status == 0)[]
  ]
})

addExperiments(
  prob.designs = list(piecewise_constant_hazard = design),
  algo.designs = list(models_fit = data.table(
    cut_step = 0.1, alpha = 0.05, maxdepth = 2, min_events = 20,
    K_max = 1, goodman_levels = "goodman"
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[n == 1000 & delta == 1][1, job.id], external = TRUE)
ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 50)]
submitJobs(ids)
waitForJobs()
