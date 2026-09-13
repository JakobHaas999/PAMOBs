### Framework for simulating data that comes from a PEM with covariates

sim_pem_data <- function(n,
                         covariate_spec = NULL,
                         formula = NULL,
                         lambdas,
                         cuts,
                         beta = NULL,
                         censoring_rate = NULL) {
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
