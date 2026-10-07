source("setup.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("simulation/scenarios/01_level_shift.R")

# Parameters ------------------------------------------------------------------

n_validation <- 10000
tau <- 5
admin_time <- 10
interval <- 0.05
censoring_rate <- 0.4
covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -1, max = 1),
  x3 = list(distfun = rbinom, size = 1, prob = 0.5)
)
covariate_effects <- c(x1 = 0.3, x2 = 0.2, x3 = 0.5)

validation_design <- expand.grid(
  delta = c(0, 0.8),
  covariate_setting = c("none", "included"),
  stringsAsFactors = FALSE
)
validation_design$seed <- 9000 + seq_len(nrow(validation_design))

beta0 <- -2

simulate <- function(delta, spec, effects) {
  scenario_level_shift(
    n = n_validation, delta = delta, tau = tau, beta0 = beta0,
    censoring_rate = censoring_rate, covariate_spec = spec,
    covariate_effects = effects, admin_time = admin_time,
    sim_interval = interval
  )
}

baseline_truth <- c("(Intercept)" = beta0)
true_hazard <- function(t, delta) exp(beta0 + delta * I(t >= tau))
cut <- sort(unique(c(seq(0, admin_time, by = interval), tau, admin_time)))

# Validation ------------------------------------------------------------------

validation_01 <- rbindlist(lapply(seq_len(nrow(validation_design)), function(i) {
  cfg <- validation_design[i, , drop = FALSE]
  included <- cfg$covariate_setting == "included"
  spec <- if (included) covariate_spec else NULL
  effects <- if (included) covariate_effects else NULL
  covariates <- if (included) names(effects) else character()

  set.seed(cfg$seed)
  dat <- simulate(cfg$delta, spec, effects)

  ped <- pammtools::as_ped(
    reformulate(covariates, response = "Surv(time, status)"),
    data = dat,
    cut = cut
  )

  ped$post_tau <- as.integer(ped$tstart >= tau)

  fit <- glm(
    reformulate(
      c("post_tau", covariates, "offset(offset)"),
      response = "ped_status"
    ),
    family = poisson(link = "log"),
    data = ped
  )

  truth <- c(baseline_truth, post_tau = cfg$delta, effects)
  estimates <- coef(summary(fit))
  stopifnot(all(names(truth) %in% rownames(estimates)))

  starts <- head(cut, -1)
  newdata <- data.frame(
    tstart = starts,
    post_tau = as.integer(starts >= tau),
    offset = 0
  )
  for (variable in covariates) newdata[[variable]] <- 0

  hazard_hat <- predict(fit, newdata = newdata, type = "response")
  hazard_true <- true_hazard(starts, cfg$delta)
  hazard_ise <- sum(diff(cut) * (hazard_hat - hazard_true)^2)
  newdata$hazard_hat <- hazard_hat
  newdata$hazard_true <- hazard_true

  data.table(
    delta = cfg$delta,
    covariate_setting = cfg$covariate_setting,
    truth = truth,
    estimate = estimates[names(truth), "Estimate"],
    se = estimates[names(truth), "Std. Error"],
    hazard_ise = hazard_ise,
    event_rate = mean(dat$status == 1),
    random_censoring = mean(dat$censoring_reason == "random")
  )
}))

# Results ---------------------------------------------------------------------
print(validation_01)
