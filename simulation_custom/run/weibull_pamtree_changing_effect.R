########################################################
## Pamtree with time-varying effects ------------------
#######################################################

### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pamtree.R")
library("parallel")

n_sim <- 200
n_obs <- 1000
n_pilot <- 10000

beta1 <- 0.4
beta2 <- 0.6

deltas <- c(0, 0.3, 0.6, 0.9)
tau <- 4

admin_time <- 10
censoring_rate <- 0.4
alpha <- 0.05

n_cores <- 2

log_baseline_hazard <- function(t) log(1.5) + log(0.05) + 0.5 * log(t)

# The formula environment holds delta, so each scenario keeps its own value
make_sim_formula <- function(delta) {
  ~ log_baseline_hazard(t) + beta1 * x1 + beta2 * x2 + delta * x1 * (t >= tau)
}
sim_formulas <- lapply(deltas, make_sim_formula)

gen_data <- function(n) {
  data.frame(
    id = seq_len(n),
    x1 = rnorm(n),
    x2 = runif(n, min = -1, max = 1)
  )
}

fit_formula <- Surv(time, status) ~ s(tend, bs = "ps") + x1 + x2 + offset(offset) | tend
fit_cut <- seq(0, admin_time, by = 0.1)

## Calibration per delta: censoring rate and smoothing parameter
set.seed(2024)
lambda_c_vals <- numeric(length(deltas))
sp_vals <- numeric(length(deltas))

for (d in seq_along(deltas)) {
  message("Calibration for delta = ", deltas[d])

  event_times_pilot <- sim_event_times(
    formula = sim_formulas[[d]],
    data = gen_data(n_pilot),
    breaks = tau,
    maxt = admin_time
  )
  lambda_c_vals[d] <- calibrate_censoring(event_times_pilot, censoring_rate, admin_time)

  sim_data_gs <- sim_surv(
    formula = sim_formulas[[d]],
    data = gen_data(n_obs),
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c_vals[d]
  )
  gs_res <- pam_grid_search(
    formula = fit_formula,
    data = sim_data_gs,
    cut = fit_cut,
    sp_space = c(0, 0.1, 1, 10, 50, 100, 1000),
    method = "cv",
    folds = 3
  )
  sp_vals[d] <- gs_res$selected_sp
}

## Simulation
design <- expand.grid(rep = seq_len(n_sim), delta_id = seq_along(deltas))
design$seed <- 8100 + design$rep
n_jobs <- nrow(design)

run_job <- function(j) {
  d <- design$delta_id[j]
  i <- design$rep[j]
  set.seed(design$seed[j])

  sim_data <- sim_surv(
    formula = sim_formulas[[d]],
    data = gen_data(n_obs),
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c_vals[d]
  )

  fit <- tryCatch(
    pamtree(
      formula = fit_formula,
      data = sim_data,
      cut = fit_cut,
      sp = sp_vals[d],
      maxdepth = 2,
      alpha = alpha
    ),
    error = function(e) e
  )
  success <- !inherits(fit, "error")
  p <- NA_real_
  split <- NA
  tau_hat <- NA_real_
  x1_before <- NA_real_
  x1_after <- NA_real_

  if (success) {
    root <- node_party(fit$tree)
    root_p <- info_node(root)$p.value
    if (length(root_p) == 1 && is.finite(root_p)) {
      p <- unname(root_p)
    }
    split <- length(nodeids(fit$tree, terminal = TRUE)) > 1
    if (split) {
      tau_hat <- breaks_split(split_node(root))
      # Rows of coef() are the terminal nodes: first = tend <= tau_hat, second = tend > tau_hat
      cf <- coef(fit$tree)
      x1_before <- cf[1, "x1"]
      x1_after <- cf[2, "x1"]
    }
  }

  file.create(file.path(progress_dir, j))
  n_done <- length(list.files(progress_dir))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf(
    "[%s] delta = %.1f, rep %d done | %d/%d (%.0f%%) | %.1f min elapsed, ~%.1f min left",
    format(Sys.time(), "%H:%M:%S"), deltas[d], i, n_done, n_jobs, 100 * n_done / n_jobs,
    elapsed, elapsed / n_done * (n_jobs - n_done)
  ))

  data.frame(
    rep = i,
    seed = design$seed[j],
    delta = deltas[d],
    sp = sp_vals[d],
    success = success,
    error = if (success) NA_character_ else conditionMessage(fit),
    root_p_value = p,
    split_detected = split,
    tau_hat = tau_hat,
    x1_before = x1_before,
    x1_after = x1_after,
    observed_censoring = mean(sim_data$status == 0)
  )
}

progress_dir <- file.path(tempdir(), "progress")
dir.create(progress_dir, showWarnings = FALSE)
start_time <- Sys.time()

results <- mclapply(
  seq_len(n_jobs),
  run_job,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

failed <- vapply(results, inherits, logical(1), what = "try-error")
if (any(failed)) {
  warning("Worker crashed in jobs: ", paste(which(failed), collapse = ", "))
}
results <- do.call(rbind, results[!failed])
rownames(results) <- NULL

## Save
results_dir <- "simulation_custom/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    results = results,
    calibration = data.frame(delta = deltas, lambda_c = lambda_c_vals, sp = sp_vals),
    settings = list(
      n_sim = n_sim, n_obs = n_obs, beta1 = beta1, beta2 = beta2, tau = tau,
      censoring_rate = censoring_rate, alpha = alpha, fit_cut = fit_cut,
      fit_formula = fit_formula
    )
  ),
  file.path(results_dir, "pamtree_time_changing_effect.rds")
)
