generate_covariables <- function(n, split_type) {
  x <- rnorm(n)
  if (split_type == "binary") {
    g <- rbinom(n, size = 1, prob = 0.5)
    z <- factor(g, levels = c(0, 1))
  } else if (split_type == "numeric") {
    z <- round(runif(n, min = -1, max = 1), 2)
    g <- as.integer(z > 0)
  } else {
    stop("Unknown split_type")
  }

  data.frame(
    id = seq_len(n),
    x = x,
    z = z,
    g = g
  )
}

generate_data <- function(n, scenario, split_type, delta, beta, admin_time, lambda_c) {
  data <- generate_covariables(n, split_type)

  log_baseline <- switch(scenario,
    cubic = function(t) -2.5 + 0.25 * t - 0.025 * t^2 + 0.001 * t^3,
    weibull = function(t) log(0.05) + log(1.5) + 0.5 * log(t),
    stop("Unknown scenario")
  )

  sim_formula <- ~ log_baseline(t) + beta * x + delta * g * x

  sim_surv(
    formula = sim_formula,
    data = data,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
}

fit_and_summarise <- function(instance, bs_df, cut_step, alpha, maxdepth, min_events) {
  t0 <- proc.time()[["elapsed"]]
  data <- instance[, c("time", "status", "x", "z")]

  fit <- tryCatch(
    pemtree(
      formula = Surv(time, status) ~ bs(tend, df = bs_df) + x + offset(offset) | z,
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
      split = NA, p_root = NA_real_, split_point = NA_real_,
      n_terminal = NA_integer_, n_events = sum(data$status),
      runtime = proc.time()[["elapsed"]] - t0,
      error = conditionMessage(fit)
    ))
  }

  root <- node_party(fit$tree)
  split <- !is.terminal(root)
  brk <- if (split) breaks_split(split_node(root)) else NULL

  list(
    split = split,
    p_root = info_node(root)$p.value,
    split_point = if (is.null(brk)) NA_real_ else brk,
    n_terminal = length(nodeids(fit$tree, terminal = TRUE)),
    n_events = sum(data$status),
    runtime = proc.time()[["elapsed"]] - t0,
    error = NA_character_
  )
}

calibrate_one <- function(scenario, split_type, delta, beta,
                          admin_time, censoring_rates, n_calibration) {
  d <- generate_data(
    n = n_calibration, scenario = scenario, split_type = split_type,
    delta = delta, beta = beta, admin_time = admin_time, lambda_c = 0
  )
  event_times <- data.frame(eventtime = d$time)

  data.table::data.table(
    censoring_rate = censoring_rates,
    lambda_c = vapply(
      censoring_rates, function(rate) {
        calibrate_censoring(event_times, censoring_rate = rate, admin_time = admin_time)
      },
      numeric(1)
    )
  )
}
