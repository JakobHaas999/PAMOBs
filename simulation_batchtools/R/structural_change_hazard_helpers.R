generate_data <- function(n, scenario, delta, tau, admin_time, lambda_c) {
  data <- data.frame(id = seq_len(n))

  log_baseline <- switch(scenario,
    cubic = function(t) -2.5 + 0.25 * t - 0.025 * t^2 + 0.001 * t^3,
    weibull = function(t) log(0.05) + log(1.5) + 0.5 * log(t),
    stop("Unknown scenario")
  )

  sim_formula <- ~ log_baseline(t) + delta * (t >= tau)

  sim_surv(
    formula = sim_formula,
    data = data,
    breaks = tau,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
}

calibrate_one <- function(n_calibration, scenario, delta, tau, admin_time, censoring_rate) {
  d <- generate_data(
    n = n_calibration,
    scenario = scenario,
    delta = delta,
    tau = tau,
    admin_time = admin_time,
    lambda_c = 0
  )

  event_times <- data.frame(eventtimes = d$time)

  data.table(
    censoring_rate = censoring_rate,
    lambda_c = vapply(censoring_rate, function(rate) {
      calibrate_censoring(event_times, rate, admin_time)
    }, numeric(1))
  )
}

fit_and_summarise <- function(instance, cut_step, alpha, maxdepth, min_events) {
  data <- instance[, c("time", "status")]

  fit_args <- list(
    data = data,
    min_events = min_events,
    alpha = alpha,
    maxdepth = maxdepth
  )

  if (!is.na(cut_step)) {
    fit_args$cut <- seq(0, ceiling(max(data$time)), by = cut_step)
  }

  fits <- lapply(c(4, 6, 8), function(bs_df) {
    tryCatch(
      do.call(pemtree, c(
        list(formula = Surv(time, status) ~ bs(tend, df = bs_df) + offset(offset) | tend),
        fit_args
      )),
      error = function(e) e
    )
  })


  summarise_fit <- function(fit, df) {
    if (inherits(fit, "error")) {
      return(data.table(
        df = df,
        split = NA,
        p_root = NA_real_,
        tau_hat = NA_real_,
        error = conditionMessage(fit)
      ))
    }

    root <- node_party(fit$tree)
    split <- !is.terminal(root)
    p <- info_node(root)$p.value

    data.table(
      df = df,
      split = split,
      p_root = if (length(p) == 1) unname(p) else NA_real_,
      tau_hat = if (split) breaks_split(split_node(root)) else NA_real_,
      error = NA_character_
    )
  }

  out <- rbindlist(Map(summarise_fit, fits, c(4, 6, 8)))
  out[, observed_censoring := mean(instance$status == 0)]
  out
}
