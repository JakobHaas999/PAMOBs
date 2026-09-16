### Simulate covariates for a PEM or PAM
#
## Example call
# set.seed(123)
# n <- 1000
# spec <- list(
#   age = list(
#     distfun = rnorm,
#     mean = 44,
#     sd = 8
#   ),
#   treatment = list(
#     distfun = rbinom,
#     size = 1,
#     prob = 0.5
#   )
# )
# x <- sim_covariates(n, spec)
# head(x)
#> age treatment
#> 1 39.51619         0
#> 2 42.15858         0
#> 3 56.46967         0
#> 4 44.56407         1
#> 5 45.03430         0
#> 6 57.72052         1
sim_covariates <- function(n, spec) {
  # Arguments:
  # n: Number of simulations/observations
  # spec: List of Specifications of covariates

  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_list(spec, min.len = 1)

  covariates <- lapply(spec, function(s) {
    checkmate::assert_list(s, names = "unique")
    checkmate::assert_true("distfun" %in% names(s))

    distfun <- s$distfun
    args <- s[names(s) != "distfun"]

    do.call(distfun, c(list(n = n), args))
  })

  as.data.frame(covariates)
}

### Simulate covariates for a PEM or PAM
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
sim_pem <- function(n, lambdas, cuts = NULL, X = NULL, beta = NULL) {
  if (length(cuts) == 1 && cuts == 0) {
    cuts <- NULL
  }

  # asserts
  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_numeric(
    lambdas,
    any.missing = FALSE,
    lower = 0,
    min.len = 1
  )
  checkmate::assert_true(all(lambdas > 0))
  checkmate::assert_numeric(cuts,
    any.missing = FALSE,
    len = length(lambdas) - 1,
    null.ok = TRUE
  )
  if (length(cuts) > 1) {
    checkmate::assert_true(all(diff(cuts) > 0))
  }
  if (!is.null(X)) {
    checkmate::assert(
      checkmate::check_numeric(X, any.missing = FALSE, len = n),
      checkmate::check_matrix(X, any.missing = FALSE, nrows = n),
      combine = "or"
    )
    n_features <- if (is.matrix(X)) ncol(X) else 1
    checkmate::assert_numeric(beta, any.missing = FALSE, len = n_features)
    X <- as.matrix(X)
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


### Simulate censoring times with a given rate
#
## Example function call
# set.seed(123)
# times <- sim_pem(
#   n = 1000,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5)
# )
# censoring_rate <- 0.6
# surv_data <- sim_censoring(times, censoring_rate)
# head(surv_data)
# #> time status
# #> 1 1.0250338      0
# #> 2 0.5394712      0
# #> 3 9.5810974      1
# #> 4 0.3157736      1
# #> 5 0.5621098      1
# #> 6 1.4354894      0
# mean(surv_data$status == 0)
# #> [1] 0.604
#
## Check beta estimates with censored data
# set.seed(234)
# n <- 1000
# x <- rbinom(n, 1, prob = 0.5)
# true_beta <- log(2)
# times <- sim_pem(
#   n = n,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5),
#   X = matrix(x, ncol = 1),
#   beta = true_beta
# )
# surv_data <- sim_censoring(times, 0.6)
# surv_data <- cbind(surv_data, x)
# ped <- pammtools::as_ped(
#   survival::Surv(time, status) ~ x,
#   data = surv_data,
#   cut = c(2, 5)
# )
# fit <- glm(
#   formula = ped_status ~ interval + x,
#   data = ped,
#   offset = ped$offset,
#   family = poisson(link = "log")
# )
# exp(coef(fit)[["x"]])
# #> [1] 2.007998

sim_censoring <- function(times, censoring_rate) {
  # Arguments:
  # times: numeric vector of event times
  # censoring_rate: Proportion of censored observations

  # asserts
  checkmate::assert_numeric(times, any.missing = FALSE, lower = 0)
  checkmate::assert_true(all(times > 0))
  checkmate::assert_number(censoring_rate)
  checkmate::assert_true(censoring_rate > 0 && censoring_rate < 1)

  lambda_c <- uniroot(
    f = function(lambda_c) {
      p_censored <- mean(1 - exp(-lambda_c * times))
      p_censored - censoring_rate
    },
    interval = c(1e-8, 10000)
  )$root

  censoring_times <- rexp(length(times), rate = lambda_c)

  data.frame(
    time = pmin(times, censoring_times),
    status = as.integer(times <= censoring_times)
  )
}


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
