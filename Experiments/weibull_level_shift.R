### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

n_sim <- n_obs <- 1000
n_pilot <- 10000

admin_time <- 10
censoring_rate <- 0.4
alpha <- 0.05

cut <- seq(0, admin_time, by = 0.1)

sim_formula <- ~ log(0.05) + log(1.5) + 0.5 * log(t) +
  delta * (t >= tau)
deltas <- log(c(1.25, 1.5, 2))
tau <- 3

# Calirate lambda censoring
lambda_c <- numeric(length(deltas))
calibration_seeds <- c(101, 102, 103)

for (i in seq_along(deltas)) {
  cat(sprintf("Running calibration: %s/%s", i, length(deltas)), "\n")
  set.seed(calibration_seeds[i])

  delta <- deltas[i]
  pilot_events <- sim_event_times(
    formula = sim_formula,
    data = data.frame(id = seq_len(n_pilot)),
    breaks = tau,
    maxt = admin_time
  )
  lambda_c[i] <- calibrate_censoring(
    event_times =pilot_events,
    censoring_rate = censoring_rate,
    admin_time = admin_time
  )
}

fit_formulas <- list(
  spline_time = Surv(time, status) ~ bs(tend, df = 6) + offset(offset) | tend,
  reference_time = Surv(time, status) ~ log(tend) + offset(offset) | tend
)

results <- vector("list", n_sim * length(deltas) * length(fit_formulas))
seeds <- matrix(
  1710 + seq_len(n_sim * length(deltas)), nrow = n_sim, ncol = length(deltas)
)

k <- 0

for (j in seq_along(deltas)) {
  delta <- deltas[j]
  current_lambda_c <- lambda_c[j]
  cat("\n\n\nRunning simulations for delta =", exp(delta), "\n\n")

  for (i in seq_len(n_sim)) {
    cat("Simulation", i, "of", n_sim, "\n")
    seed <- seeds[i, j]
    set.seed(seed)

    sim_data <- sim_surv(
      formula = sim_formula,
      data = data.frame(id = seq_len(n_obs)),
      breaks = tau,
      admin_time = admin_time,
      lambda_c = current_lambda_c
    )

    for (model_name in names(fit_formulas)) {
      k <- k + 1
      fit <- tryCatch(
        pemtree(
          formula = fit_formulas[[model_name]],
          data = sim_data,
          cut = cut,
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
      }


      results[[k]] <- data.frame(
        reep = i,
        seeds = seeds[i, j],
        model = model_name,
        success = success,
        error = if (success) NA_character_ else conditionMessage(fit),
        root_p_value = p,
        split_detected = split,
        observed_censoring = mean(sim_data$status == 0)
      )
    }
  }
}

results <- do.call(rbind, results)
rownames(results) <- NULL
# Save results
results_dir <- "Experiments/Results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(results_dir, "weibul_level_shift.rds"))