sim_event_times <- function(
  formula,
  data,
  breaks = numeric(0),
  interval = c(1e-8, 500),
  maxt = NULL,
  integration_tol = 1e-8
) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  if (length(formula) != 2L) {
    stop("formula must be one-sided.")
  }

  rhs <- formula[[2L]]
  formula_env <- environment(formula)

  loghazard <- function(t, x, betas, ...) {
    values <- as.list(x)
    values$t <- t

    eta <- eval(rhs, envir = values, enclos = formula_env)
    rep_len(eta, length.out = length(t))
  }

  args <- list(x = data, idvar = "id", ids = data$id, interval = interval, maxt = maxt)

  if (!length(breaks)) {
    args$loghazard <- loghazard
  } else {
    breaks <- sort(unique(breaks))
    args$cumhazard <- function(t, x, betas, ...) {
      vapply(t, function(end) {
        if (end <= 0) {
          return(0)
        }

        knots <- c(0, breaks[breaks > 0 & breaks < end], end)
        sum(vapply(seq_len(length(knots) - 1), function(j) {
          integrate(
            function(s) exp(loghazard(s, x, betas, ...)),
            lower = knots[j], upper = knots[j + 1],
            rel.tol = integration_tol, abs.tol = integration_tol
          )$value
        }, FUN.VALUE = numeric(1)))
      }, FUN.VALUE = numeric(1))
    }
  }

  simulated <- do.call(simsurv::simsurv, args)
  simulated[match(data$id, simulated$id), , drop = FALSE]
}

calibrate_censoring <- function(event_times, censoring_rate, admin_time = Inf) {
  checkmate::assert_number(censoring_rate, lower = 0, upper = 1)
  if (censoring_rate == 1) stop("censoring_rate has to be lower than 1")
  checkmate::assert_number(admin_time, lower = 0)

  event_times <- event_times$eventtime
  admin_censoring <- if (is.finite(admin_time)) {
    mean(event_times >= admin_time)
  } else {
    0
  }
  if (censoring_rate <= admin_censoring) {
    if (censoring_rate < admin_censoring) {
      warning(
        sprintf(
          paste0(
            "Target censoring rate (%.3f) is below the administrative",
            "censoring rate (%.3f). No additional random censoring is applied."
          ),
          censoring_rate, admin_censoring
        )
      )
    }
    return(0)
  }

  objective <- function(lambda) {
    mean(
      1 - (event_times < admin_time) * exp(-lambda * event_times)
    ) - censoring_rate
  }

  upper <- 1
  while (objective(upper) < 0) {
    upper <- upper * 2
    if (!is.finite(upper)) stop("Could not bracket censoring rate")
  }

  uniroot(
    objective,
    interval = c(0, upper), tol = 1e-10
  )$root
}

sim_censoring <- function(event_times, lambda_c, admin_time = Inf) {
  checkmate::assert_number(lambda_c, lower = 0)
  checkmate::assert_number(admin_time, lower = 0)

  n <- length(event_times)

  censor_times <- if (lambda_c == 0) {
    rep(Inf, n)
  } else {
    rexp(n, rate = lambda_c)
  }

  time <- pmin(event_times, censor_times, admin_time)
  status <- as.integer(
    event_times <= censor_times & event_times < admin_time
  )

  data.frame(time = time, status = status)
}

sim_surv <- function(
  formula,
  data,
  breaks = numeric(0),
  admin_time = Inf,
  lambda_c = 0,
  interval = NULL,
  integration_tol = 1e-8
) {
  checkmate::assert_number(admin_time, lower = 0)
  if (admin_time <= 0) stop(" admin_time must be positive.")
  
  if (!"id" %in% colnames(data)) data$id <- seq_len(nrow(data))
  
  maxt <- if (is.finite(admin_time)) admin_time else NULL
  if (is.null(interval)) {
    upper <- if (is.null(maxt)) 500 else 1.01 * maxt
    interval <- c(1e-8, upper)
  }

  events <- sim_event_times(
    formula = formula,
    data = data,
    breaks = breaks,
    interval = interval,
    maxt = maxt,
    integration_tol = integration_tol
  )
  if (
    anyNA(events$eventtime) ||
      any(!is.finite(events$eventtime)) ||
      anyNA(events$status) ||
      any(events$status != 1L)
  ) {
    stop(
      "Some event times were not successfully simulated. ",
      "Check the simulation interval."
    )
  }
  censored <- sim_censoring(event_times = events$eventtime, lambda_c = lambda_c, admin_time = admin_time)

  cbind.data.frame(data, censored)
}

## Example call
# shape <- 1.5
# rate <- 0.05
# delta <- log(3)
# tau <- 5
#
# log_hazard <- ~ log(rate) + log(shape) + (shape - 1) * log(t) + delta * (t >= tau)
# data <- data.frame(id = seq_len(3000))
#
# set.seed(3)
# sim_data_3 <- sim_surv(
#   formula = log_hazard, data = data, breaks = tau, censoring_rate = 0.3
# )
#
# truth <- data.frame(t = sort(unique(c(seq(0, 10, by = 0.01), tau))))
# truth$hazard <- with(truth, rate * shape * t^(shape - 1) * exp(delta * (t >= tau)))
#
# cuts <- sort(unique(c(
#   seq(0, ceiling(max(sim_data_3$time)), by = 1), tau
# )))
# ped_3 <- pammtools::as_ped(
#   formula = Surv(time, status) ~ 1, data = sim_data_3, cut = cuts
# )
# data.table::setDT(ped_3)
# hazard_est_3 <- ped_3[, .(
#   events = sum(ped_status), exposure = sum(exp(offset)), n_risk = .N
# ), by = .(tstart, tend)]
# hazard_est_3[, hazard := events / exposure]
#
# ggplot() +
#   geom_segment(
#     data = hazard_est_3[tstart < 10],
#     aes(x = tstart, xend = tend, y = hazard, yend = hazard),
#     colour = "steelblue", linewidth = 0.9
#   ) +
#   geom_line(data = truth, aes(x = t, y = hazard)) +
#   geom_vline(xintercept = tau, lty = 2) +
#   labs(x = "time", y = expression(lambda(t))) +
#   theme_bw()
