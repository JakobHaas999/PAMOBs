make_level_shift_input <- function(pars, n = pars$n) {
  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_number(pars$delta, finite = TRUE)
  checkmate::assert_number(pars$tau, finite = TRUE, lower = 0)
  dat <- data.frame(id = seq_len(n))

  log_h0 <- switch(
    pars$scenario,
    cubic = function(t) {
      pars$beta0 + pars$beta1 * t + pars$beta2 * t^2 + pars$beta3 * t^3
    },
    weibull = function(t) {
      log(pars$rate) + log(pars$shape) + (pars$shape - 1) * log(t)
    },
    lognormal = function(t) {
      stats::dlnorm(t, meanlog = pars$meanlog, sdlog = pars$sdlog, log = TRUE) -
        stats::plnorm(t, meanlog = pars$meanlog, sdlog = pars$sdlog, lower.tail = FALSE, log.p = TRUE)
    },
    stop("Unknown scenario: ", pars$scenario)
  )

  if (pars$covariate_setting == "included") {
    checkmate::assert_number(pars$gamma1, finite = TRUE)
    checkmate::assert_number(pars$gamm2, finite = TRUE)
    dat$x1 <- stats::rnorm(n)
    dat$x2 <- stats::rbinom(n, size = 1, prob = 0.5)

    formula <- ~ log_h0(t) + pars$delta * (t >= pars$tau) +
      pars$gamma1 * x1 + pars$gamm2 * x2
  } else if (pars$covariate_setting == "none") {
    formula <- ~ log_h0(t) + pars$delta * (t >= pars$tau)
  } else {
    stop("Unknown covariate setting: ", pars$covariate_setting)
  }

  list(data = dat, formula = formula)
}

generate_level_shift_data <- function(pars) {
  checkmate::assert_number(pars$lambda_c, finite = TRUE, lower = 0)

  input <- make_level_shift_input(pars)
  sim_surv(
    formula = input$formula,
    data = input$data,
    breaks = if (pars$delta == 0) numeric(0) else pars$tau,
    admin_time = pars$admin_time,
    lambda_c = pars$lambda_c
  )
}
