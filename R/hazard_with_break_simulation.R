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

#### Simulate from a complex hazard with a structural change
sim_rejection_break <- function(n, hazard, M) {
  # Arguments:
  # n: number of simulations

  # asserts
  checkmate::assertCount(n)
  checkmate::assertFunction(hazard, nargs = 1)
  checkmate::assertNumber(M)

  res <- numeric(n)

  for (i in seq_len(n)) {
    t <- 0

    repeat {
      w <- exp(M)
      t <- t + w
      lambda_t <- hazard(t)

      if (runif(1) <= lambda_t / M) {
        res[[i]] <- t
      }
    }
  }

  res
}
