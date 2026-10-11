### setup
source("setup.R")
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pamob.R")

log_hazard <- function(t) {
  -1.1 * sin(pi / 4) * t
}

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = 0, max = 3),
  x3 = list(distfun = rbinom, prob = 0.5, size = 1)
)

admin_time <- 10
cut <- seq(0, admin_time, by = 0.05)

design <- expand.grid(
  n = c(1000, 2000),
  delta = c(0, 0.4, 0.8),
  censoring_rate = 0.3,
  rep = seq_len(10),
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 8000 + seq_len(nrow(design))

results <- lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, , drop = FALSE]
  set.seed(d$seed)

  sim_formula <- as.formula(
    paste0(
      "~ log_hazard(t) + 0.7 * x1 + 0.4 * x2 + ",
      d$delta,
      " * (x3 == 1) * x1"
    )
  )

  out <- tryCatch(
    {
      sim_data <- sim_pam(
        n = d$n,
        formula = sim_formula,
        covariate_spec = covariate_spec,
        censoring_rate = d$censoring_rate,
        admin_time = admin_time,
        cut = cut
      )

      fit <- pamob(
        formula = Surv(time, status) ~ s(tend, bs = "ps", k = 8) + x1 + x2 + offset(offset) | x3,
        data = sim_data,
        cut = cut,
        min_events = 20,
        alpha = 0.05,
        maxdepth = 2,
        vcov = "opg"
      )
      n_terminal <- length(partykit::nodeids(fit$tree, terminal = TRUE))
      root_info <- partykit::info_node(partykit::node_party(fit$tree))
      test_failed <- is.null(root_info$test) ||
        anyNA(root_info$test["p.value", ])
      if (test_failed) stop("MOB fluctuation test failed")

      data.frame(
        success = TRUE,
        split_detected = n_terminal > 1,
        n_terminal_nodes = n_terminal,
        observed_censoring = mean(sim_data$status == 0),
        error = NA_character_
      )
    },
    error = function(e) {
      data.frame(
        success = FALSE,
        split_detected = NA,
        n_terminal_nodes = NA_integer_,
        observed_censoring = NA_real_,
        error = conditionMessage(e)
      )
    }
  )

  cbind(d, out)
})

results <- do.call(rbind, results)
saveRDS(results, "simulation/results/pilot/06_covar_instability_pamob_pilot.rds")

# Inspect results
data.table::setDT(results)

# Null scenario
results[success == TRUE & delta == 0, .(
  alpha_hat = mean(split_detected),
  mcse = sqrt(mean(split_detected) * (1 - mean(split_detected)) / .N)
), by = .(n)]
# Power
results[success == TRUE & delta != 0, .(
  power_hat = mean(split_detected),
  mcse = sqrt(mean(split_detected) * (1 - mean(split_detected)) / .N)
), by = .(n, delta)]
