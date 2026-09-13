### Framework for simulating survival data from a PEM with covariates
#
## Example function call without covariates
# set.seed(123)
# surv_data <- sim_pem_data(
#   n = 1000,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5),
#   censoring_rate = 0.4
# )
# head(surv_data)
#
## Example function call with covariates
# set.seed(123)
#
# covariate_spec <- list(
#   x1 = list(
#     dist = "normal",
#     mean = 0,
#     sd = 1
#   ),
#   x2 = list(
#     dist = "binomial",
#     size = 1,
#     prob = 0.5
#   )
# )
#
# surv_data <- sim_pem_data(
#   n = 1000,
#   covariate_spec = covariate_spec,
#   formula = ~ x1 + x2,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5),
#   beta = c(0.5, -0.3),
#   censoring_rate = 0.4
# )
# head(surv_data)
#
## Example function call with an interaction
# set.seed(123)
#
# surv_data <- sim_pem_data(
#   n = 1000,
#   covariate_spec = covariate_spec,
#   formula = ~ x1 * x2,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5),
#   beta = c(
#     x1 = 0.5,
#     x2 = -0.3,
#     `x1:x2` = 0.8
#   ),
#   censoring_rate = 0.4
# )
# head(surv_data)
sim_pem_data <- function(n,
                         covariate_spec = NULL,
                         formula = NULL,
                         lambdas,
                         cuts,
                         beta = NULL,
                         censoring_rate = NULL) {
  # Arguments:
  # n: number of individuals to simulate
  # covariate_spec: specification passed to sim_covariates();
  #   NULL if no covariates should be simulated
  # formula: one-sided formula defining the covariate effects,
  #   e.g. ~ x1 + x2 or ~ x1 * x2
  # lambdas: piecewise constant baseline hazard rates
  # cuts: cut points defining the time intervals of the baseline hazard
  # beta: regression coefficients corresponding to the columns of the
  #   model matrix generated from formula
  # censoring_rate: desired proportion of censored observations;
  #   NULL means no censoring

  if (is.null(covariate_spec) != is.null(formula)) {
    stop("'covariate_spec' and 'formula' must either both be NULL or both be supplied.")
  }
  if (!is.null(formula) && is.null(beta)) {
    stop("'beta' must be supplied when 'formula' is supplied.")
  }
  if (is.null(formula) && !is.null(beta)) {
    stop("'beta' must be NULL when no 'formula' is supplied.")
  }

  # Simulate covariates
  if (!is.null(covariate_spec)) {
    covariates <- sim_covariates(n = n, spec = covariate_spec)
    X <- model.matrix(formula, covariates)
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  } else {
    covariates <- NULL
    X <- NULL
  }

  # Check beta
  if (!is.null(X)) {
    if (length(beta) != ncol(X)) {
      stop(
        sprintf(
          "beta has length %s, but formula produces %s model matrix columns",
          length(beta), ncol(X)
        )
      )
    }
  }

  # Simulate event times
  times <- sim_pem(
    n = n,
    lambdas = lambdas,
    cuts = cuts,
    X = X, beta = beta
  )

  # Simulate censoring times
  if (!is.null(censoring_rate)) {
    surv_data <- sim_censoring(times = times, censoring_rate = censoring_rate)
  } else {
    surv_data <- data.frame(
      time = times,
      status = rep(1L, n)
    )
  }

  cbind.data.frame(surv_data, covariates)
}
