### Simulating from a PEM
#
## Example validation without covariates:
# set.seed(123)
# t <- sim_pem(
#   n = 100000,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5)
# )
# mean(t > 2)
# exp(-0.2)
# mean(t > 5)
# exp(-(0.1 * 2 + 0.3 * 3))
#
## Example with covariates
# set.seed(234)
# n <- 10000
# X <- cbind(
#   x1 = rnorm(n, mean = 2, sd = 1),
#   x2 =  rbinom(n, 1, 0.4)
# )
# beta <- c(0.5, -0.2)
# t <- sim_pem(
#   n = n,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5),
#   X = X,
#   beta = beta
# )
# head(t)
sim_pem <- function(n, lambdas, cuts, X = NULL, beta = NULL) {
  # TODO comments

  # asserts
  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_numeric(
    lambdas,
    any.missing = FALSE,
    lower = 0,
    min.len = 1
  )
  checkmate::assert_true(all(lambdas > 0))
  checkmate::assert_numeric(cuts, any.missing = FALSE, len = length(lambdas) - 1)
  if (length(cuts) > 1) {
    checkmate::assert_true(all(diff(cuts) > 0))
  }
  if (!is.null(X)) {
    checkmate::assert_matrix(X, any.missing = FALSE, nrows = n)
    checkmate::assert_numeric(beta, any.missing = FALSE, len = ncol(X))
  }

  cuts <- c(0, cuts, Inf)
  times <- numeric(n)

  for (i in seq_len(n)) {
    # linear predictor
    if (is.null(X)) {
      eta_i <- 0
    } else {
      eta_i <- sum(X[i, ] * beta)
    }
    lambdas_i <- lambdas * exp(eta_i)

    E <- rexp(1)
    H_before <- 0

    for (j in seq_along(lambdas)) {
      int_length <- cuts[[j + 1]] - cuts[[j]]
      H_next <- H_before + int_length * lambdas_i[[j]]

      if (E <= H_next) {
        t <- (E - H_before) / lambdas_i[[j]] + cuts[[j]]
        times[[i]] <- t
        break
      }
      H_before <- H_next
    }
  }

  times
}
