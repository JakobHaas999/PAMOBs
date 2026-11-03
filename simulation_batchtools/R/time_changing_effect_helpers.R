generate_data <- function(n, scenario, delta, tau, beta, admin_time, lambda_c) {
  data <- data.frame(id = seq_len(n), x = rnorm(n))

  log_baseline <- switch(scenario,
    cubic = function(t) -2.5 + 0.25 * t - 0.025 * t^2 + 0.001 * t^3,
    weibull = function(t) log(0.05) + log(1.5) + 0.5 * log(t),
    stop("Unknown scenario")
  )

  sim_formula <- ~ log_baseline(t) + beta * x + delta * x * (t >= tau)

  sim_surv(
    formula = sim_formula,
    data = data,
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
}

calibrate_one <- function(scenario, delta, tau, beta, admin_time, censoring_rate, n_calibration) {
  d <- generate_data(
    n = n_calibration, scenario = scenario, delta = delta, tau = tau, beta = beta,
    admin_time = admin_time, lambda_c = 0
  )
  event_times <- data.frame(eventtime = d$time)

  data.table::data.table(
    censoring_rate = censoring_rate,
    lambda_c = vapply(
      censoring_rate, function(rate) {
        calibrate_censoring(event_times, rate, admin_time)
      }, numeric(1)
    )
  )
}


fit_and_summarise <- function(instance, bs_df, cut_step, alpha, maxdepth, min_events) {
  t0 <- proc.time()[["elapsed"]]
  data <- instance[, c("time", "status", "x")]

  fit <- tryCatch(
    pemtree(
      formula = Surv(time, status) ~ bs(tend, df = bs_df) + x + offset(offset) | tend,
      data = data,
      cut = seq(0, ceiling(max(data$time)), by = cut_step),
      min_events = min_events,
      alpha = alpha,
      maxdepth = maxdepth
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      split = NA,
      p_root = NA_real_,
      split_point = NA_real_,
      root_valid = NA,
      n_terminal = NA_integer_,
      x_before = NA_real_,
      x_after = NA_real_,
      n_events = sum(data$status),
      runtime = proc.time()[["elapsed"]] - t0,
      error = conditionMessage(fit)
    ))
  }

  root <- node_party(fit$tree)
  split <- !is.terminal(root)
  brk <- if (split) breaks_split(split_node(root)) else NULL
  p <- info_node(root)$p.value
  cf <- coef(fit$tree)

  list(
    split = split,
    p_root = if (length(p) == 1) unname(p) else NA_real_,
    split_point = if (is.null(brk)) NA_real_ else brk,
    root_valid = is.finite(info_node(root)$objfun),
    n_terminal = length(nodeids(fit$tree, terminal = TRUE)),
    x_before = if (split) cf[1, "x"] else NA_real_,
    x_after = if (split) cf[2, "x"] else NA_real_,
    n_events = sum(data$status),
    runtime = proc.time()[["elapsed"]] - t0,
    error = NA_character_
  )
}
