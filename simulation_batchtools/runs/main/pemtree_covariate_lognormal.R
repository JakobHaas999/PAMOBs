Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
### Comparison of Weibull AFT tree and PEM tree with a lognormal baseline
library("batchtools")
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

# ----- 1. Design ------
run_id <- "pemtree_covariate_lognormal_v2"
n_rep <- 1000
beta <- 0.6
admin_time <- 10
n_calibration <- 50000

design <- CJ(
  n = c(250, 500, 1000),
  delta = seq(0, 1, by = 0.25),
  censoring_rate = c(0.3, 0.5, 0.7)
)
design[, `:=`(beta = beta, admin_time = admin_time)]

# ----- 2. Calibration -----
set.seed(67)
calib <- unique(design[, .(delta, censoring_rate)])[,
  {
    df <- data.frame(
      id = seq_len(n_calibration),
      x = rnorm(n_calibration),
      z = round(runif(n_calibration, -1, 1), 2)
    )
    lognormal <- function(t) {
      dlnorm(t, meanlog = log(4), sdlog = 1.25, log = TRUE) -
        plnorm(t, meanlog = log(4), sdlog = 1.25, lower.tail = FALSE, log.p = TRUE)
    }
    d <- sim_surv(
      formula = ~ lognormal(t) + beta * x + delta * x * (z > 0),
      data = df,
      admin_time = admin_time,
      lambda_c = 0
    )
    event_times <- data.frame(eventtime = d$time)
    data.table(
      censoring_rate = censoring_rate,
      lambda_c = vapply(censoring_rate, function(r) {
        calibrate_censoring(event_times, r, admin_time = admin_time)
      }, numeric(1))
    )
  },
  by = "delta"
]
n_before <- nrow(design)
design <- merge(design, calib, by = c("delta", "censoring_rate"))
stopifnot(nrow(design) == n_before)

# ----- 3. Registry -----
reg_dir <- "simulation_batchtools/registries/main"
dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, run_id),
  seed = 6757,
  packages = c(
    "survival", "partykit", "pammtools", "splines", "simsurv",
    "Formula", "checkmate", "sandwich", "data.table"
  ),
  source = c("R/sim_surv.R", "R/pemtree_splines.R")
)
reg$cluster.functions <- makeClusterFunctionsMulticore(ncpus = 16)

addProblem("covariate", fun = function(data, job, n, delta, beta,
                                       admin_time, lambda_c, ...) {
  df <- data.frame(
    id = seq_len(n),
    x = rnorm(n, mean = 0, sd = 1),
    z = round(runif(n, min = -1, max = 1), 2)
  )
  lognormal <- function(t) {
    dlnorm(t, meanlog = log(4), sdlog = 1.25, log = TRUE) -
      plnorm(t, meanlog = log(4), sdlog = 1.25, lower.tail = FALSE, log.p = TRUE)
  }

  sim_surv(
    formula = ~ lognormal(t) + beta * x + delta * x * (z > 0),
    data = df,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
})

addAlgorithm("fit_models", fun = function(data, job, instance, cut_step, alpha,
                                          maxdepth, min_events, ...) {
  data <- instance[, c("time", "status", "x", "z")]
  fits <- vector("list", 2)

  make_wbreg_fit <- function(min_events) {
    force(min_events)
    function(y, x, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
      invalid <- list(
        coefficients = setNames(numeric(ncol(x)), colnames(x)),
        objfun = Inf, estfun = NULL, object = NULL
      )
      if (sum(y[, "status"]) < min_events) {
        return(invalid)
      }

      fit <- tryCatch(
        suppressWarnings(survreg(y ~ 0 + x, weights = weights, dist = "weibull")),
        error = function(e) NULL
      )

      if (is.null(fit) || any(!is.finite(coef(fit)))) {
        return(invalid)
      }

      list(
        coefficients = coef(fit),
        objfun = -as.numeric(logLik(fit)),
        estfun = if (estfun) estfun(fit) else NULL,
        object = if (object) fit else NULL
      )
    }
  }

  fits[[1]] <- tryCatch(
    pemtree(
      formula = Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | z,
      data = data,
      cut = seq(0, ceiling(max(data$time)), by = cut_step),
      min_events = min_events,
      alpha = alpha,
      maxdepth = maxdepth
    ),
    error = function(e) e
  )

  fits[[2]] <- tryCatch(
    mob(
      formula = Surv(time, status) ~ x | z,
      data = data,
      fit = make_wbreg_fit(min_events = min_events),
      control = mob_control(alpha = alpha, maxdepth = maxdepth)
    ),
    error = function(e) e
  )

  summarise_fit <- function(fit, tree, model) {
    if (inherits(fit, "error")) {
      return(data.table(
        model = model, split = NA, p_root = NA_real_, split_point = NA_real_,
        error = conditionMessage(fit)
      ))
    }

    root <- node_party(tree)
    split <- !is.terminal(root)
    p <- info_node(root)$p.value
    brk <- if (split) breaks_split(split_node(root)) else NULL
    data.table(
      model = model, split = split,
      p_root = if (length(p) == 1) p else NA_real_,
      split_point = if (is.null(brk)) NA_real_ else brk,
      error = NA_character_
    )
  }

  rbind(
    summarise_fit(fits[[1]], if (!inherits(fits[[1]], "error")) fits[[1]]$tree, "pemtree"),
    summarise_fit(fits[[2]], fits[[2]], "weibull_aft")
  )[, observed_censoring := mean(instance$status == 0)][]
})

addExperiments(
  prob.designs = list(covariate = design),
  algo.designs = list(fit_models = data.table(
    cut_step = 0.1, alpha = 0.05,
    maxdepth = 2, min_events = 20
  )),
  repls = n_rep
)

pars <- unwrap(getJobPars())
testJob(pars[n == 1000 & delta == 0.75][1, job.id], external = TRUE)

ids <- findNotDone()
ids[, chunk := chunk(job.id, chunk.size = 25)]
submitJobs(ids)
waitForJobs()
