###############################################################################
## Significance level check pemtree_splines ###################################
###############################################################################

### setup

# packages
source("setup.R")
# simulation functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
# pemtree class functions
source("R/pemtree_splines.R")
source("R/pemtree.R")

# --- Scenario 1: Constant baseline hazard function --------------------------------
# ----------------------------------------------------------------------------------

# Specify parameters
bazehaz <- function(t, beta0) beta0
n <- 1000
beta0 <- -2
covariate_spec <- list(
  x1 = list(
    distfun = rnorm,
    mean = 2,
    sd = 2
  ),
  x2 = list(
    distfun = runif,
    min = -2,
    max = 2
  ),
  x3 = list(
    distfun = rbinom,
    prob = 0.5,
    size = 1
  )
)
formula <- ~ bazehaz(t, beta0 = beta0) + 0.2 * x1 + 0.08 * x2 + 0.6 * x3
censoring_rate <- 0.4
admin_time <- 10
cut <- seq(0, max(admin_time), by = 0.05)

n_sim <- 1000
seeds <- 222 + seq_len(n_sim)

results <- data.frame(
  seed = seeds,
  success = FALSE,
  split = NA_integer_,
  error = NA_character_
)

cat("Null simulation for picewise constant baseline hazard:\n\n\n")
for (i in seq_len(n_sim)) {
  cat("Running repetition: ", i, "\n\n")
  set.seed(seeds[i])

  result <- tryCatch(
    {
      sim_data <- sim_pam(
        n = n,
        formula = formula,
        covariate_spec = covariate_spec,
        censoring_rate = censoring_rate,
        admin_time = admin_time,
        cut = cut
      )

      fit <- pemtree_splines(
        formula = Surv(time, status) ~ x1 + x2 + x3 + offset(offset) | tend,
        data = sim_data,
        cut = cut,
        min_events = 30,
        maxdepth = 2,
        alpha = 0.05
      )

      list(
        success = TRUE,
        split = as.integer(
          length(nodeids(fit$tree, terminal = TRUE)) > 1L
        ),
        error = NA_character_
      )
    },
    error = function(e) {
      list(
        success = FALSE,
        split = NA_integer_,
        error = conditionMessage(e)
      )
    }
  )

  results$success[i] <- result$success
  results$split[i] <- result$split
  results$error[i] <- result$error
}

alpha_hat <- mean(results$split[results$success])
mcse <- sqrt(alpha_hat * (1 - alpha_hat) / sum(results$success))
test <- binom.test(
  sum(results$split[results$success]),
  sum(results$success),
  p = 0.05
)
ci_lower <- test$conf.int[1]
ci_upper <- test$conf.int[2]
c(
  alpha_hat = alpha_hat,
  mcse = mcse,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)


# --- Scenario 2:  Flexible baseline hazard function -------------------------------
# ----------------------------------------------------------------------------------

# Specify parameters
basehaz <- function(t) {
  -2.5 + 0.25 * t - 0.025 * t^2
}
n <- 1000
covariate_spec <- list(
  x1 = list(
    distfun = rnorm,
    mean = 2,
    sd = 2
  ),
  x2 = list(
    distfun = runif,
    min = -2,
    max = 2
  ),
  x3 = list(
    distfun = rbinom,
    prob = 0.5,
    size = 1
  )
)
formula <- ~ basehaz(t) + 0.2 * x1 + 0.08 * x2 + 0.6 * x3
censoring_rate <- 0.4
admin_time <- 10
cut <- seq(0, admin_time, by = 0.05)

n_sim <- 1000
seeds <- 22222 + seq_len(n_sim)

results_splines <- data.frame(
  seed = seeds,
  success = FALSE,
  split = NA_integer_,
  error = NA_character_
)

cat("Null simulation for flexible hazard:\n\n\n")
for (i in seq_len(n_sim)) {
  cat("Running repetition: ", i, "\n\n")
  set.seed(results_splines$seed[i])

  result <- tryCatch(
    {
      sim_data <- sim_pam(
        n = n,
        formula = formula,
        covariate_spec = covariate_spec,
        censoring_rate = censoring_rate,
        admin_time = admin_time,
        cut = cut
      )

      fit <- pemtree_splines(
        formula = Surv(time, status) ~ bs(tend, df = 5, degree = 3) +
          x1 + x2 + x3 + offset(offset) | tend,
        data = sim_data,
        cut = seq(0, admin_time, by = 0.1),
        min_events = 50,
        alpha = 0.05,
        maxdepth = 2
      )

      list(
        success = TRUE,
        split = as.integer(
          length(nodeids(fit$tree, terminal = TRUE)) > 1L
        ),
        error = NA_character_
      )
    },
    error = function(e) {
      list(
        success = FALSE,
        split = NA_integer_,
        error = conditionMessage(e)
      )
    }
  )
  results_splines$success[i] <- result$success
  results_splines$split[i] <- result$split
  results_splines$error[i] <- result$error
}

alpha_hat_splines <- mean(results_splines$split[results_splines$success])
mcse_splines <- sqrt(alpha_hat_splines * (1 - alpha_hat_splines) / sum(results_splines$success))
test_splines <- binom.test(
  x = sum(results_splines$split[results_splines$success]),
  n = sum(results_splines$success),
  p = 0.05
)
ci_lower_splines <- test_splines$conf.int[1]
ci_upper_splines <- test_splines$conf.int[2]
c(
  alpha_hat = alpha_hat_splines,
  mcse = mcse_splines,
  ci_lower = ci_lower_splines,
  ci_upper = ci_upper_splines
)
