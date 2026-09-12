### Framework for simulating data that comes from a PEM with covariates

sim_pem_data <- function(n,
                         covariate_spec = NULL,
                         lambdas,
                         cuts,
                         beta = NULL,
                         censoring_rate = NULL) {
  if (is.null(covariate_spec) != is.null(beta)) {
    stop("covariate_spec and beta must either both be NULL or both be supplied.")
  }

  if (!is.null(covariate_spec)) {
    X <- sim_covariates(n = n, spec = covariate_spec)
  } else {
    X <- NULL
  }

  times <- sim_pem(
    n = n,
    lambdas = lambdas,
    cuts = cuts,
    X = X, beta = beta
  )

  if (!is.null(censoring_rate)) {
    surv_data <- sim_censoring(times = times, censoring_rate = censoring_rate)
  } else {
    surv_data <- data.frame(
      time = times,
      status = rep(1L, n)
    )
  }

  cbind.data.frame(surv_data, X)
}
