source("setup.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("simulation/scenarios/02_level_shift_splines.R")

# Parameters ------------------------------------------------------------------

n_validation <- 10000
tau <- 5
admin_time <- 10
interval <- 0.05
censoring_rate <- 0.4

basehaz <- function(t) {
  -2.5 + 0.25 * t - 0.025 * t^2
}

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
validation_design$seed <- 10000 + seq_len(nrow(validation_design))

simulate <- function(delta, spec, effects) {
  scenario_level_shift_splines(
    n = n_validation,
    delta = delta,
    tau = tau,
    basehaz_fun = basehaz,
    censoring_rate = censoring_rate,
    covariate_spec = spec,
    covariate_effects = effects,
    admin_time = admin_time,
    sim_interval = interval
  )
}

true_hazard <- function(t, delta) {
  exp(basehaz(t) + delta * as.integer(t >= tau))
}

baseline_truth <- c(
  "(Intercept)" = -2.5,
  "tstart" = 0.25,
  "I(tstart^2)" = -0.025
)

cut <- sort(unique(c(
  seq(0, admin_time, by = interval),
  tau,
  admin_time
)))
starts <- head(cut, -1)

# Validation ------------------------------------------------------------------

validation_runs <- lapply(seq_len(nrow(validation_design)), function(i) {
  cfg <- validation_design[i, , drop = FALSE]
  included <- cfg$covariate_setting == "included"
  spec <- if (included) covariate_spec else NULL
  effects <- if (included) covariate_effects else NULL
  covariates <- if (included) names(effects) else character()

  set.seed(cfg$seed)
  dat <- simulate(cfg$delta, spec, effects)

  stopifnot(
    nrow(dat) == n_validation,
    !anyNA(dat),
    all(dat$time > 0 & dat$time <= admin_time),
    all(dat$status %in% c(0L, 1L))
  )

  ped <- pammtools::as_ped(
    reformulate(covariates, response = "Surv(time, status)"),
    data = dat,
    cut = cut
  )
  ped$post_tau <- as.integer(ped$tstart >= tau)

  fit <- glm(
    reformulate(
      c(
        "tstart", "I(tstart^2)", "post_tau",
        covariates, "offset(offset)"
      ),
      response = "ped_status"
    ),
    family = poisson(link = "log"),
    data = ped
  )

  if (
    !isTRUE(fit$converged) ||
      fit$rank < length(coef(fit)) ||
      any(!is.finite(coef(fit)))
  ) {
    stop("Invalid oracle fit in validation case ", i)
  }

  truth <- c(baseline_truth, post_tau = cfg$delta, effects)
  estimates <- coef(summary(fit))

  newdata <- data.frame(
    tstart = starts,
    post_tau = as.integer(starts >= tau),
    offset = 0
  )
  for (variable in covariates) newdata[[variable]] <- 0

  newdata$hazard_hat <- predict(fit, newdata, type = "response")
  newdata$hazard_true <- true_hazard(starts, cfg$delta)

  hazard_ise <- sum(
    diff(cut) * (newdata$hazard_hat - newdata$hazard_true)^2
  )

  plot <- ggplot(newdata, aes(x = tstart)) +
    geom_step(
      aes(y = hazard_true, linetype = "Truth"),
      direction = "hv"
    ) +
    geom_step(
      aes(y = hazard_hat, linetype = "Estimate"),
      direction = "hv"
    ) +
    labs(
      title = paste(
        "delta =", cfg$delta,
        "| covariates =", cfg$covariate_setting
      ),
      subtitle = sprintf("Hazard ISE = %.6f", hazard_ise),
      x = "t",
      y = expression(lambda(t)),
      linetype = NULL
    ) +
    theme_bw()

  list(
    summary = data.table(
      delta = cfg$delta,
      covariate_setting = cfg$covariate_setting,
      term = names(truth),
      truth = unname(truth),
      estimate = estimates[names(truth), "Estimate"],
      se = estimates[names(truth), "Std. Error"],
      hazard_ise = hazard_ise,
      event_rate = mean(dat$status == 1L),
      random_censoring = mean(dat$censoring_reason == "random"),
      administrative_censoring =
        mean(dat$censoring_reason == "administrative")
    ),
    plot = plot
  )
})

# Results ---------------------------------------------------------------------

validation_02 <- rbindlist(
  lapply(validation_runs, `[[`, "summary")
)

print(validation_02)

for (result in validation_runs) {
  print(result$plot)
}
