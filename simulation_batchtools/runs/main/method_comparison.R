### Comparison of pemtree and exposure-scaled pemtree
### setup
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/pemtree_exact.R")

# ----- 1. Design -----
run_id <- "method_comparison_v1"
n_rep <- 1000
beta0 <- -2
tau <- 3
admin_time <- 10
n_calibration <- 50000

design <- CJ(
  n = c(250, 500),
  delta = seq(0, 1, by = 0.25),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(beta0 = beta0, tau = tau, admin_time = admin_time)]

# ----- 2. Calibration -----
set.seed(302342)
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

# ----- 3. Registry -----
reg_path <- "simulation_batchtools/registries/main"
dir.create(reg_path, showWarnings = FALSE, recursive = TRUE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_path, run_id),
  seed = 32223,
  packages = c(
    "survival", "partykit", "pammtools", "simsurv", "Formula", "checkmate",
    "sandwich", "data.table"
  ),
  source = c(
    "R/sim_surv.R", "R/pemtree_splines.R", "R/pemtree_exact.R"
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

addAlgorithm("models_fit", fun = function(data, job, instance, alpha, maxdepth,
                                          min_events, pi0, R, ...) {
  data <- instance[, c("time", "status")]

  fits <- vector("list", length = 4)

  pemtree_args <- list(
    formula = Surv(time, status) ~ offset(offset) | tend,
    data = data,
    min_events = min_events,
    alpha = alpha,
    maxdepth = maxdepth,
    trim = pi0
  )

  fits[[1]] <- tryCatch(
    {
      pemtree_args$cut <- seq(0, ceiling(max(data$time)), by = 0.1)
      do.call(pemtree, pemtree_args)
    },
    error = function(e) e
  )
  fits[[2]] <- tryCatch(
    {
      pemtree_args$cut <- NULL
      do.call(pemtree, pemtree_args)
    },
    error = function(e) e
  )
  fits[[3]] <- tryCatch(
    {
      ped <- as_ped(Surv(time, status) ~ 1, data = data, cut = seq(0, ceiling(max(data$time)), by = 0.1))
      pem_tree(ped, alpha = alpha, pi0 = pi0, R = R, minevents = min_events, maxdepth = maxdepth)
    },
    error = function(e) e
  )
  fits[[4]] <- tryCatch(
    {
      ped <- as_ped(Surv(time, status) ~ 1, data = data)
      pem_tree(ped, alpha = alpha, pi0 = pi0, R = R, minevents = min_events, maxdepth = maxdepth)
    },
    error = function(e) e
  )

  summarise_pemtree <- function(fit, label) {
    if (inherits(fit, "error")) {
      return(data.table(
        model = label, split = NA, p_root = NA_real_, split_point = NA_real_,
        lambda1 = NA_real_, lambda2 = NA_real_, error = conditionMessage(fit)
      ))
    }
    root <- node_party(fit$tree)
    split <- !is.terminal(root)
    p <- info_node(root)$p.value
    data.table(
      model = label, split = split,
      p_root = if (length(p) == 1) unname(p) else NA_real_,
      split_point = if (split) breaks_split(split_node(root)) else NA_real_,
      lambda1 = if (split) exp(coef(fit$tree)[["2"]]) else NA_real_,
      lambda2 = if (split) exp(coef(fit$tree)[["3"]]) else NA_real_,
      error = NA_character_
    )
  }

  summarise_pem_tree <- function(fit, label) {
    if (inherits(fit, "error")) {
      return(data.table(
        model = label, split = NA, p_root = NA_real_, split_point = NA_real_,
        lambda1 = NA_real_, lambda2 = NA_real_, error = conditionMessage(fit)
      ))
    }
    split <- !fit[1, ]$terminal
    p <- fit[1, ]$p.value
    data.table(
      model = label, split = split,
      p_root = p,
      split_point = if (split) fit[1, ]$split else NA_real_,
      lambda1 = if (split) fit[2, ]$lambda else NA_real_,
      lambda2 = if (split) fit[3, ]$lambda else NA_real_,
      error = NA_character_
    )
  }

  labels <- c("pemtree_grid", "pemtree_event", "es_grid", "es_event")
  out <- rbindlist(lapply(seq_along(fits), function(i) {
    if (i <= 2) summarise_pemtree(fits[[i]], labels[i]) else summarise_pem_tree(fits[[i]], labels[i])
  }))
  out[, observed_censoring := mean(data$status == 0)]
  out[]
})

addExperiments(
  prob.designs = list(piecewise_constant_hazard = design),
  algo.designs = list(models_fit = data.table(
    alpha = 0.05, maxdepth = 2,
    min_events = 20, pi0 = 0.1, R = 1e4
  )),
  repls = n_rep
)

# ----- 4. Test and submit jobs -----
pars <- unwrap(getJobPars())
testJob(pars[n == 500 & delta == 1 & censoring_rate == 0.3][1, job.id], external = TRUE)
ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 50)]
submitJobs(ids)
waitForJobs()
