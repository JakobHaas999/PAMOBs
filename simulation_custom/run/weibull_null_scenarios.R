### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")

n_sim <- 1000
n_obs <- 1000
n_pilot <- 10000

admin_time <- 10
censoring_rate <- 0.4
alpha <- 0.05

cut <- seq(0, admin_time, by = 0.1)

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

fit_formulas <- list(
  spline_time = Surv(time, status) ~ bs(tend, df = 6) + offset(offset) | tend,
  reference_time = Surv(time, status) ~ log(tend) + offset(offset) | tend,
  spline_noise = Surv(time, status) ~ bs(tend, df = 6) + offset(offset) | z_noise,
  reference_noise = Surv(time, status) ~ log(tend) + offset(offset) | z_noise
)

seeds <- 5100 + seq_len(n_sim)
results <- vector("list", n_sim * length(fit_formulas))
k <- 0

for (i in seq_len(n_sim)) {
  cat("\nSimulation", i, "of", n_sim, "\n")
  set.seed(seeds[i])

  sim_data <- sim_surv(
    formula = sim_formula,
    data = data.frame(id = seq_len(n_obs)),
    admin_time = admin_time,
    lambda_c = lambda_c
  )
  sim_data$z_noise <- runif(n_obs)

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

    if (success) {
      root <- node_party(fit$tree)
      root_p <- info_node(root)$p.value

      if (length(root_p) == 1 && is.finite(root_p)) {
        p <- root_p
      }
      split <- length(nodeids(fit$tree, terminal = TRUE)) > 1
    }

    results[[k]] <- data.frame(
      rep = i,
      seed = seeds[[i]],
      model = model_name,
      success = success,
      error = if (success) NA_character_ else conditionMessage(fit),
      root_p_value = p,
      split_detected = split,
      observed_censoring = mean(sim_data$status == 0)
    )
  }
}

results <- do.call(rbind, results)
rownames(results) <- NULL
## Save
results_dir <- "simulation_custom/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(results_dir, "weibull_null_scenarios.rds"))
