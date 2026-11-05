generate_experiment <- function(setting, n, delta, tau, beta, lambda_c, admin_time) {
  log_baseline <- function(t) log(0.05) + log(1.5) + 0.5 * log(t)

  if (setting == "covariate") {
    x <- rnorm(n)
    z <- round(runif(n, min = -1, max = 1), 2)
    g <- as.integer(z > 0)

    data <- data.frame(id = seq_len(n), x = x, z = z, g = g)
    sim_formula <- ~ log_baseline(t) + beta * x + delta * g * x

    return(sim_surv(
      formula = sim_formula,
      data = data,
      admin_time = admin_time,
      lambda_c = lambda_c
    ))
  } else if (setting == "time_changing_effect") {
    x <- rnorm(n)
    data <- data.frame(id = seq_len(n), x = x)

    sim_formula <- ~ log_baseline(t) + beta * x + delta * x * (t >= tau)

    return(sim_surv(
      formula = sim_formula,
      data = data,
      breaks = tau,
      admin_time = admin_time,
      lambda_c = lambda_c
    ))
  } else if (setting == "structural_change") {
    data <- data.frame(id = seq_len(n))

    sim_formula <- ~ log_baseline(t) + delta * (t >= tau)

    return(sim_surv(
      formula = sim_formula,
      data = data,
      breaks = tau,
      admin_time = admin_time,
      lambda_c = lambda_c
    ))
  } else {
    stop("Unknown setting")
  }
}

calibrate_one <- function(setting, n_calibration, delta, tau, beta, admin_time, censoring_rate) {
  d <- generate_experiment(
    setting = setting,
    n = n_calibration,
    delta = delta,
    tau = tau,
    beta = beta,
    lambda_c = 0,
    admin_time = admin_time
  )
  event_times <- data.frame(eventtimes = d$time)

  data.table::data.table(
    censoring_rate = censoring_rate,
    lambda_c = vapply(censoring_rate, function(rate) {
      calibrate_censoring(event_times, rate, admin_time)
    }, numeric(1))
  )
}

fit_and_summarise <- function(instance, setting, cut_step, admin_time, alpha, min_events, B, frac) {
  formula <- switch(setting,
    covariate = Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | z,
    time_changing_effect = Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | tend,
    structural_change = Surv(time, status) ~ bs(tend, df = 4) + offset(offset) | tend,
    stop("Unknown setting")
  )

  vars <- switch(setting,
    covariate = c("time", "status", "x", "z"),
    time_changing_effect = c("time", "status", "x"),
    structural_change = c("time", "status")
  )
  data <- instance[, vars]

  fit_args <- list(
    formula = formula,
    data = data,
    cut = seq(0, admin_time, by = cut_step),
    alpha = alpha,
    min_events = min_events
  )

  res <- split_stability(pemtree, fit_args = fit_args, B = B, frac = frac)

  tree_full <- res$full_fit$tree
  root <- node_party(tree_full)
  full_split <- !is.terminal(root)
  full_split_point <- if (full_split) breaks_split(split_node(root)) else NA_real_
  full_split_var <- if (full_split) colnames(tree_full$data)[split_node(root)$varid] else NA_character_
  p <- info_node(root)$p.value
  full_p_value <- if (length(p) == 1) p else NA_real_

  data.table(
    b = seq_len(B),
    split = res$split,
    split_var = res$split_var,
    split_point = res$split_point,
    p_value = res$p_value,
    full_split = rep(full_split, B),
    full_split_var = rep(full_split_var, B),
    full_split_point = rep(full_split_point, B),
    full_p_value = rep(full_p_value, B),
    observed_censoring = rep(mean(instance$status == 0), B)
  )
}

# # Test
# for (s in c("covariate", "time_changing_effect", "structural_change")) {
#   d <- generate_experiment(s, n = 1000, delta = 1, tau = 3, beta = 0.6, lambda_c = 0.05, admin_time = 10)
#   print(fit_and_summarise(
#     d, s,
#     cut_step = 0.1, admin_time = 10, alpha = 0.05, min_events = 20, B = 5, frac = 0.632
#   ))
# }
