### Function to generate event times from a Weibull distribution
# set.seed(12)
# shape <- 2
# rate <- 0.02
# delta <- log(3)
# tau <- 5
#
# dat <- sim_weibull(
#   n = 10000,
#   shape = shape,
#   rate = rate,
#   censoring_rate = 0.3,
#   tau = tau,
#   delta = delta
# )
#
# ped <- pammtools::as_ped(survival::Surv(time, status) ~ 1, data = dat, cut = seq(0, 10, by = 0.25))
# data.table::setDT(ped)
# hazard_est <- ped[, .(
#   tend = tend[1],
#   n_risk = .N,
#   events = sum(ped_status),
#   exposure = sum(exp(offset)),
#   hazard = sum(ped_status) / sum(exp(offset))
# ), .(interval)]
#
# plot(hazard_est$tend, hazard_est$hazard,
#   type = "s",
#   xlab = "t", ylab = expression(hat(lambda(t)))
# )
# abline(v = tau, lty = 2)


sim_weibull <- function(
  n,
  shape,
  rate,
  covariate_spec = NULL,
  covariate_effects = NULL,
  censoring_rate,
  tau = NULL,
  delta = 0
) {
  # Arguments:
  # n: number of simulations
  # shape: shape parameter of weibull distribution
  # rate: rate parameter of weibull distribution
  # covariate_spec:
  # covariate_effects:
  # censoring_rate:
  # tau:
  # delta:

  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_number(shape, lower = .Machine$double.eps)
  checkmate::assert_number(rate, lower = .Machine$double.eps)
  checkmate::assert_number(censoring_rate, lower = 0, upper = 1)
  if (!is.null(tau)) {
    checkmate::assert_number(tau, lower = 0)
    checkmate::assert_number(delta)
  }

  if (is.null(covariate_spec)) {
    if (!is.null(covariate_effects)) {
      warning("covariate_effects should be NULL if covariate_spec is NULL")
    }

    eta <- rep(0, n)
    X <- NULL
  } else {
    checkmate::assert_list(covariate_spec, names = "unique")

    if (is.null(covariate_effects)) {
      message("covariate_effects set to 1")
      covariate_effects <- rep(1, length(covariate_spec))
      names(covariate_effects) <- names(covariate_spec)
    }

    checkmate::assert_numeric(
      covariate_effects,
      len = length(covariate_spec)
    )

    checkmate::assert_names(
      names(covariate_effects),
      permutation.of = names(covariate_spec)
    )

    covariate_effects <- covariate_effects[names(covariate_spec)]

    data_list <- lapply(covariate_spec, function(spec) {
      checkmate::assert_true("distfun" %in% names(spec))

      fun <- spec$distfun
      args <- spec[names(spec) != "distfun"]

      do.call(fun, c(list(n = n), args))
    })

    X <- do.call(cbind, data_list)
    colnames(X) <- names(covariate_spec)

    eta <- drop(X %*% covariate_effects)
  }

  # Generate event times
  e <- -log(runif(n))

  if (is.null(tau)) {
    t <- (
      e / (rate * exp(eta))
    )^(1 / shape)
  } else {
    H_tau <- rate * tau^shape * exp(eta)


    before_tau <- e <= H_tau
    t <- numeric(n)

    t[before_tau] <- (
      e[before_tau] /
        (rate * exp(eta[before_tau]))
    )^(1 / shape)

    t[!before_tau] <- (
      tau^shape +
        (e[!before_tau] - H_tau[!before_tau]) /
          (rate * exp(eta[!before_tau] + delta))
    )^(1 / shape)
  }

  # Censoring rate
  if (censoring_rate == 0) {
    ct <- rep(Inf, n)
  } else {
    objective <- function(lambda_c) {
      mean(1 - exp(-lambda_c * t)) - censoring_rate
    }

    lambda_c <- uniroot(
      f = objective,
      interval = c(1e-10, 1e4)
    )$root

    ct <- rexp(n, rate = lambda_c)
  }
  status <- as.integer(t <= ct)
  time <- pmin(t, ct)

  if (!is.null(covariate_spec)) {
    out <- data.frame(time = time, status = status, X)
  } else {
    out <- data.frame(time = time, status = status)
  }
  out$id <- seq_len(nrow(out))
  return(out)
}
