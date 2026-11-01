########################################################
## Pamtree demo: effect of mob_control(parm = ...) ----
#######################################################
# Single runs (no repetitions) on two DGPs with a break at tau:
#   baseline_break: baseline hazard jumps at tau, covariate effects constant
#   effect_break:   effect of x1 changes at tau, baseline smooth
# Each data set is fitted with pamtree once testing all parameters
# (parm = NULL) and once testing only the covariate effects.

### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pamtree.R")

n_obs <- 1000
n_pilot <- 5000
tau <- 4
admin_time <- 10
censoring_rate <- 0.4
alpha <- 0.05
sp <- 10

beta1 <- 0.4
beta2 <- 0.6
delta_baseline <- log(2)
delta_x1 <- 0.9

log_baseline_hazard <- function(t) log(1.5) + log(0.05) + 0.5 * log(t)

sim_formulas <- list(
  baseline_break = ~ log_baseline_hazard(t) + beta1 * x1 + beta2 * x2 +
    delta_baseline * (t >= tau),
  effect_break = ~ log_baseline_hazard(t) + beta1 * x1 + beta2 * x2 +
    delta_x1 * x1 * (t >= tau)
)

parms <- list(
  all = NULL,
  covariates = c("x1", "x2")
)

fit_formula <- Surv(time, status) ~ s(tend, bs = "ps") + x1 + x2 + offset(offset) | tend
fit_cut <- seq(0, admin_time, by = 0.1)

gen_data <- function(n) {
  data.frame(
    id = seq_len(n),
    x1 = rnorm(n),
    x2 = runif(n, min = -1, max = 1)
  )
}

results <- list()

for (dgp in names(sim_formulas)) {
  set.seed(4711)
  pilot <- sim_event_times(
    formula = sim_formulas[[dgp]],
    data = gen_data(n_pilot),
    breaks = tau,
    maxt = admin_time
  )
  lambda_c <- calibrate_censoring(pilot, censoring_rate, admin_time)

  sim_data <- sim_surv(
    formula = sim_formulas[[dgp]],
    data = gen_data(n_obs),
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )

  for (parm_name in names(parms)) {
    message("DGP: ", dgp, " | parm: ", parm_name)
    start <- Sys.time()
    fit <- pamtree(
      formula = fit_formula,
      data = sim_data,
      cut = fit_cut,
      sp = sp,
      maxdepth = 2,
      alpha = alpha,
      parm = parms[[parm_name]]
    )
    runtime <- as.numeric(difftime(Sys.time(), start, units = "secs"))

    root <- node_party(fit$tree)
    split <- !is.null(split_node(root))
    cf <- coef(fit$tree)

    results[[length(results) + 1]] <- data.frame(
      dgp = dgp,
      parm = parm_name,
      root_p_value = unname(info_node(root)$p.value),
      split_detected = split,
      tau_hat = if (split) breaks_split(split_node(root)) else NA_real_,
      x1_before = if (split) cf[1, "x1"] else NA_real_,
      x1_after = if (split) cf[2, "x1"] else NA_real_,
      runtime_sec = runtime
    )
    print(fit$tree)
  }
}

results <- do.call(rbind, results)
print(results, digits = 3)
