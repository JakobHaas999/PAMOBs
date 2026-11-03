generate_data <- function(n, admin_time, lambda_c) {
  data <- data.frame(id = seq_len(n))

  log_hazard <- function(t) log(1.5) + log(0.05) + 0.5 * log(t)

  sim_surv(
    formula = ~ log_hazard(t),
    data = data,
    admin_time = admin_time,
    lambda_c = lambda_c
  )
}

fit_and_summarise <- function(instance, cut_step, alpha, maxdepth, min_events) {
  data <- instance[, c("time", "status")]

  fits <- vector("list", 4)
  # names(fits) <-
  fit_args <- list(
    data = data,
    cut = seq(0, ceiling(max(data$time)), by = cut_step),
    min_events = min_events,
    alpha = alpha,
    maxdepth = maxdepth
  )
  fits[[1]] <- tryCatch(
    {
      formula <- Surv(time, status) ~ log(tend) + offset(offset) | tend
      do.call(pemtree, c(list(formula = formula), fit_args))
    },
    error = function(e) e
  )
  fits[[2]] <- tryCatch(
    {
      formula <- Surv(time, status) ~ s(tend, k = 10, bs = "ps") + offset(offset) | tend
      sp <- 0
      do.call(pamtree, c(list(formula = formula, sp = sp), fit_args))
    },
    error = function(e) e
  )
  fits[[3]] <- tryCatch(
    {
      formula <- Surv(time, status) ~ s(tend, k = 10, bs = "ps") + offset(offset) | tend
      sp <- 30
      do.call(pamtree, c(list(formula = formula, sp = sp), fit_args))
    },
    error = function(e) e
  )
  fits[[4]] <- tryCatch(
    {
      formula <- Surv(time, status) ~ s(tend, bs = "ps", k = 10) + offset(offset) | tend
      gam_method <- "REML"
      do.call(pamtree, c(list(formula = formula, gam_method = gam_method), fit_args))
    },
    error = function(e) e
  )

  summarise_fit <- function(fit, variant) {
    if (inherits(fit, "error")) {
      return(data.table(
        variant = variant, split = NA, p_root = NA_real_, tau_hat = NA_real_,
        sp_root = NA_real_, error = conditionMessage(fit)
      ))
    }

    root <- node_party(fit$tree)
    split <- !is.terminal(root)
    data.table(
      variant = variant, split = split, p_root = info_node(root)$p.value,
      tau_hat = if (split) breaks_split(split_node(root)) else NA_real_,
      sp_root = if (variant == "REML") refit.modelparty(fit$tree, node = 1)$sp else NA_real_,
      error = NA_character_
    )
  }

  rbindlist(Map(summarise_fit, fits, c("oracle_log", "sp0", "sp30", "REML")))
}

calibrate_one <- function(n_calibration, admin_time, censoring_rate) {
  d <- generate_data(n_calibration, admin_time, lambda_c = 0)
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
