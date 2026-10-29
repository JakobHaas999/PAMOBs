### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
library(parallel)

## Checks whether the inflated type I error of reference_time in the null
## scenario comes from evaluating log(t) at the interval end:
##   reference_time      log(tend),  cut width 0.1  (as in weibull_null_scenarios.R)
##   reference_mid_time  log(tmid),  cut width 0.1
##   reference_time_cut  log(tend),  cut width 0.5
##   reference_time_all  log(tend),  all unique event times as cut points

n_sim <- 200
n_obs <- 1000
n_pilot <- 10000

admin_time <- 10
censoring_rate <- 0.4
alpha <- 0.05

n_cores <- 6

sim_formula <- ~ log(0.05) + log(1.5) + 0.5 * log(t)

set.seed(101)

pilot_events <- sim_event_times(
  formula = sim_formula,
  data = data.frame(id = seq_len(n_pilot)),
  maxt = admin_time
)
lambda_c <- calibrate_censoring(
  event_times = pilot_events,
  censoring_rate = censoring_rate,
  admin_time = admin_time
)

fit_specs <- list(
  reference_time = list(
    formula = Surv(time, status) ~ log(tend) + offset(offset) | tend,
    cut = seq(0, admin_time, by = 0.1)
  ),
  reference_mid_time = list(
    formula = Surv(time, status) ~ log((tstart + tend) / 2) + offset(offset) | tend,
    cut = seq(0, admin_time, by = 0.1)
  ),
  reference_time_cut = list(
    formula = Surv(time, status) ~ log(tend) + offset(offset) | tend,
    cut = seq(0, admin_time, by = 0.5)
  ),
  # cut = NULL: as_ped uses all unique event times as cut points
  reference_time_all = list(
    formula = Surv(time, status) ~ log(tend) + offset(offset) | tend,
    cut = NULL
  )
)

# Same seeds as weibull_null_scenarios.R, so the data sets match its first n_sim runs
seeds <- 5100 + seq_len(n_sim)

run_rep <- function(i) {
  message("Simulation ", i, " of ", n_sim, " started")
  set.seed(seeds[i])

  sim_data <- sim_surv(
    formula = sim_formula,
    data = data.frame(id = seq_len(n_obs)),
    admin_time = admin_time,
    lambda_c = lambda_c
  )

  rep_results <- lapply(names(fit_specs), function(model_name) {
    fit <- tryCatch(
      pemtree(
        formula = fit_specs[[model_name]]$formula,
        data = sim_data,
        cut = fit_specs[[model_name]]$cut,
        maxdepth = 2,
        alpha = alpha
      ),
      error = function(e) e
    )
    success <- !inherits(fit, "error")
    p <- NA_real_
    split <- NA
    tau_hat <- NA_real_

    if (success) {
      root <- node_party(fit$tree)
      root_p <- info_node(root)$p.value

      if (length(root_p) == 1 && is.finite(root_p)) {
        p <- unname(root_p)
      }
      split <- length(nodeids(fit$tree, terminal = TRUE)) > 1
      if (split) {
        tau_hat <- breaks_split(split_node(root))
      }
    }

    data.frame(
      rep = i,
      seed = seeds[[i]],
      model = model_name,
      success = success,
      error = if (success) NA_character_ else conditionMessage(fit),
      root_p_value = p,
      split_detected = split,
      tau_hat = tau_hat,
      observed_censoring = mean(sim_data$status == 0)
    )
  })

  # Workers do not share memory, so finished reps are counted via marker files
  file.create(file.path(progress_dir, i))
  n_done <- length(list.files(progress_dir))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  message(sprintf(
    "[%s] Simulation %d done | %d/%d (%.0f%%) | %.1f min elapsed, ~%.1f min left",
    format(Sys.time(), "%H:%M:%S"), i, n_done, n_sim, 100 * n_done / n_sim,
    elapsed, elapsed / n_done * (n_sim - n_done)
  ))

  do.call(rbind, rep_results)
}

progress_dir <- file.path(tempdir(), "progress")
dir.create(progress_dir, showWarnings = FALSE)
start_time <- Sys.time()

# Seeds are set inside run_rep, so results do not depend on n_cores.
# mc.preschedule = FALSE balances the load, since fits with a split take longer.
results <- mclapply(
  seq_len(n_sim),
  run_rep,
  mc.cores = n_cores,
  mc.preschedule = FALSE
)

failed <- vapply(results, inherits, logical(1), what = "try-error")
if (any(failed)) {
  warning("Worker crashed in reps: ", paste(which(failed), collapse = ", "))
}

results <- do.call(rbind, results[!failed])
rownames(results) <- NULL
## Save
results_dir <- "simulation_custom/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(results_dir, "weibull_null_midpoint.rds"))
