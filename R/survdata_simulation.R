#### Simulate from Cox model
# Example function call
# set.seed(123)
# n <- 1000
# X <- cbind(
#   x1 = rnorm(n, mean = 1, sd = 2),
#   x2 = sample(c(0, 1), n, replace = TRUE)
# )
# beta <- c(1.4, 0.6)
# lambda0 <- 0.3
# censoring_rate <- 0.6
# surv_data <- sim_coxph(n, X, beta, lambda0, censoring_rate)
# head(surv_data)
# mean(surv_data$status == 0)
# cox_fit <- survival::coxph(
#   formula = survival::Surv(time, status) ~ x1 + x2,
#   data = surv_data
# )
# summary(cox_fit)

sim_coxph <- function(n, X, beta, lambda0, censoring_rate) {
  # Arguments:
  # n: number of simulations
  # X: covariates:
  # beta: vector of coefficients
  # lambda0: constant baseline hazard
  # censoring_rate: number of censored observations

  # asserts
  checkmate::assertInt(n)
  checkmate::assertMatrix(X, any.missing = FALSE, mode = "numeric")
  checkmate::assertNumeric(beta, any.missing = FALSE, len = ncol(X))
  checkmate::assertNumber(lambda0)
  checkmate::assertNumber(censoring_rate)

  E <- rexp(n, rate = 1)
  event_rate <- lambda0 * exp(X %*% beta)
  T <- E / event_rate

  # Find censoring hazard
  lambda_c <- uniroot(
    f = function(lambda_c) {
      mean(lambda_c / (lambda_c + event_rate)) - censoring_rate
    },
    interval = c(1e-8, 100)
  )$root

  C <- rexp(n, rate = lambda_c)
  time <- pmin(T, C)
  status <- as.integer(T <= C)

  data.frame(
    time, status, X
  )
}

### Simulate from constant hazard with structural change
# Example function call
# set.seed(234)
# n <- 1000
# lambdas <- c(0.1, 0.7)
# tau <- 4
# censoring_rate <- 0.6
# surv_data <- sim_break(n, lambdas, tau, censoring_rate)
# hazard_fit <- survival::survfit(
#   formula = survival::Surv(time, status) ~ 1,
#   data = surv_data
# )
# plot(hazard_fit, fun = "cumhaz")
sim_break <- function(n, lambdas, tau, censoring_rate) {
  # arguments:
  # n: number of simulations
  # lambda: numeric vector of length 2 with the two constant hazard rates
  # tau: time of the structural change
  # censoring_rate: number of censored observations

  # asserts
  checkmate::assertInt(n)
  checkmate::assertNumeric(lambdas, any.missing = FALSE, len = 2)
  checkmate::assertNumber(tau)
  checkmate::assertNumber(censoring_rate, lower = 0, upper = 1)

  E <- rexp(n, rate = 1)

  lambda1 <- lambdas[[1]]
  lambda2 <- lambdas[[2]]

  T <- ifelse(E <= lambda1 * tau,
    E / lambda1,
    (E - lambda1 * tau) / lambda2 + tau
  )

  # Compute censoring rate
  lambda_c <- uniroot(
    f = function(lambda_c) {
      mean(1 - exp(-lambda_c * T)) - censoring_rate
    },
    interval = c(1e-8, 100)
  )$root

  C <- rexp(n, rate = lambda_c)
  time <- pmin(T, C)
  status <- as.integer(T <= C)

  data.frame(time, status)
}

#### Simulate from PAMM with interactions
# Example function call
# set.seed(345)
# n <- 1000
# beta <- c(0.3, 0.5, -0.6)
# interactions <- list(as.formula("~I(x1*x2)"), as.formula("~I(x1*x3)"))
# beta_interaction <- c(0.5, -0.8)
# cut <- seq(0, 4, by = 0.1)
# baseline_fun <- function(t) -2 + 0.6 * sin(t)
# surv_data <- sim_pamm_interaction(
#   n = n,
#   beta = beta,
#   interactions = interactions,
#   beta_interaction = beta_interaction,
#   cut = cut,
#   baseline_fun = baseline_fun
# )
# ped_data <- pammtools::as_ped(
#   formula = survival::Surv(time, status) ~ .,
#   data = surv_data
# )
# pamm_fit <- mgcv::gam(
#   formula = ped_status ~ s(tend) + x1 + x2 + x3 + x1:x2 + x1:x3,
#   data = ped_data,
#   offset = ped_data$offset,
#   family = poisson(),
#   method = "REML"
# )
# summary(pamm_fit)
sim_pamm_interaction <- function(n,
                                 beta,
                                 interactions,
                                 beta_interaction,
                                 cut,
                                 baseline_fun) {
  # arguments
  # n: number of simulations
  # beta: vector of coefficients
  # interactions: list of formula objects specifying the interactions
  # cut: intervall cut points
  # baseline_fun: log baseline hazard

  # asserts
  checkmate::assertCount(n)
  checkmate::assertNumeric(beta, any.missing = FALSE)
  checkmate::assertList(interactions, types = "formula", max.len = choose(length(beta), 2))
  checkmate::assertNumeric(beta_interaction, any.missing = FALSE, len = length(interactions))
  checkmate::assertNumeric(cut, any.missing = FALSE, lower = 0)
  checkmate::assertFunction(baseline_fun, nargs = 1)

  # Generate design matrix
  X <- matrix(NA_real_, nrow = n, ncol = length(beta))
  if (length(beta) >= 1) X[, 1] <- rnorm(n, mean = 4, sd = 2)
  if (length(beta) >= 2) X[, 2] <- rbinom(n, 1, prob = 0.6)
  if (length(beta) >= 3) {
    for (i in seq_along(beta)[-(1:2)]) {
      X[, i] <- runif(n, min = 0, max = 6)
    }
  }
  colnames(X) <- paste0("x", seq_along(beta))

  # Generate formula
  # Main effects
  main_parts <- paste0(beta, " * ", colnames(X))
  # interaction names
  interaction_names <- vapply(interactions, function(form) {
    attr(terms(form), "term.labels")
  }, character(1))
  # Interaction effects
  interaction_parts <- paste0(beta_interaction, " * ", interaction_names)

  all_parts <- c(
    main_parts,
    interaction_parts
  )

  # Generate full sim formula
  formula <- as.formula(
    paste0(
      "~baseline_fun(t) + ",
      paste0(all_parts, collapse = " + ")
    )
  )
  print(formula)
  # Generate sim_data
  sim_data <- pammtools::sim_pexp(
    formula = formula,
    data = as.data.frame(X),
    cut = cut
  )

  sim_data
}
