sim_surv <- function(
  formula,
  data,
  breaks = numeric(0),
  admin_time = Inf,
  censoring_rate = 0,
  interval = c(1e-8, 500),
  integration_tol = 1e-8
) {
  # asserts
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  if (length(formula) != 2) stop("formula has to be one-sided: ~ log_hazard(t).")

  checkmate::assert_number(censoring_rate, lower = 0, upper = 1)
  if (censoring_rate == 1) stop("censoring_rate has to be smaller than one.")

  out <- as.data.frame(data)
  if (!("id" %in% colnames(out))) out$id <- seq_len(nrow(out))

  rhs <- formula[[2]]
  formula_env <- environment(formula)

  log_hazard <- function(t, x, betas, ...) {
    values <- as.list(x)
    values$t <- t

    eta <- eval(rhs, envir = values, enclos = formula_env)
    rep_len(eta, length(t))
  }

  args <- list(
    x = out, idvar = "id", ids = out$id, interval = interval,
    maxt = if (censoring_rate == 0 && is.finite(admin_time)) admin_time else NULL
  )

  if (length(breaks)) {
    checkmate::assert_numeric(breaks, any.missing = FALSE, lower = 0, finite = TRUE)
    args$cumhazard <- function(t, x, betas, ...) {
      vapply(t, function(end) {
        if (end == 0) {
          return(0)
        }
        knots <- sort(unique(c(0, breaks[breaks > 0 & breaks < end], end)))
        sum(vapply(seq_len(length(knots) - 1), function(j) {
          integrate(
            function(s) exp(log_hazard(s, x, betas, ...)),
            lower = knots[j], upper = knots[j + 1], subdivisions = 1000L,
            rel.tol = integration_tol, abs.tol = integration_tol
          )$value
        }, numeric(1)))
      }, numeric(1))
    }
  } else {
    args$loghazard <- log_hazard
  }

  simulated <- do.call(simsurv::simsurv, args)
  simulated <- simulated[match(out$id, simulated$id), , drop = FALSE]

  censor_time <- rep(Inf, nrow(out))
  if (censoring_rate > 0) {
    objective <- function(lambda) {
      mean(-expm1(-lambda * simulated$eventtime)) - censoring_rate
    }
    upper <- 1
    while (objective(upper) < 0) {
      upper <- 2 * upper
      if (!is.finite(upper)) stop("Could not bracket the censoring rate.")
    }

    lambda_c <- uniroot(
      objective,
      interval = c(0, upper), tol = .Machine$double.eps, check.conv = TRUE
    )$root

    censor_time <- rexp(nrow(out), rate = lambda_c)
  }

  out$time <- pmin(simulated$eventtime, censor_time, admin_time)
  out$status <- as.integer(
    simulated$status == 1L & simulated$eventtime <= censor_time & simulated$eventtime <= admin_time
  )

  out
}

# # Example call:
#
# set.seed(12)
# covariates <- data.frame(
#   x1 = runif(1000, -1, 1),
#   x2 = runif(1000, -1, 1),
#   z = rbinom(1000, size = 1, prob = 0.5)
# )
# shape <- 2
# rate <- 0.02
# tau <- 5
#
# # The formula is evaluated as arithmetic, just as in sim_pexponential.
# # x1 * x2 means a product; it does not add main effects automatically.
# log_hazard <- ~ log(rate) + log(shape) + (shape - 1) * log(t) +
#   log(3) * (t >= tau) + # global hazard jump at tau
#   0.4 * x1 + # linear effect
#   0.6 * x2^2 + # nonlinear effect
#   0.3 * x1 * x2 + # interaction
#   (log(2) - 0.5 * log1p(t)) * z + # smoothly time-varying effect
#   0.5 * x1 * (t >= tau) # change in the x1 effect at tau
#
# dat <- sim_surv(
#   formula = log_hazard,
#   data = covariates,
#   breaks = tau,
#   admin_time = 10,
#   censoring_rate = 0.3
# )
# head(dat)
# mean(dat$status == 0)
#
# # Global jumps: list all discontinuities in breaks, e.g. breaks = c(3, 5).
# # Subject-specific change points are also possible:
# # covariates$tau_i <- runif(nrow(covariates), 3, 6)
# # dat <- sim_surv(
# #   ~ log(0.1) + 0.4 * x1 + log(3) * (t >= tau_i),
# #   data = covariates, breaks = ~ tau_i, admin_time = 10, seed = 12
# # )
# #
# # A time-varying covariate can be defined by a prespecified trajectory:
# #   ~ log(0.1) + 0.4 * (x1 + slope * t)
# # Draw slope once as a column of covariates. Do not draw random values inside
# # the formula. For analysis, represent the same trajectory in the PED data.
# #
# # For a spline-defined truth, define the spline function once outside the formula
# # (fixed knots and boundary knots), then call it inside the formula. mgcv::s()
# # describes a model to fit and is not a numeric simulation effect.
# #
# # censoring_rate targets random censoring before administrative censoring.
# # The observed total censoring fraction can therefore exceed censoring_rate.
# # interval bounds the numerical search, whereas admin_time ends follow-up.
