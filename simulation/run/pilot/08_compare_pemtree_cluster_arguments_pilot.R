### setup
# Packages
source("setup.R")
# Functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemtree_splines.R")

##### Pilot simulation

### Parameters
log_hazard <- function(t, delta) {
  -2 + delta * (t >= 5)
}

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -1, max = 1)
)

design <- expand.grid(
  n = c(1000, 2000),
  delta = c(0, 0.5),
  censoring_rate = 0.3,
  covariate_setting = c("none", "include"),
  rep = seq_len(20),
  admin_time = 10,
  sim_interval = 0.05,
  ped_interval = 0.1,
  KEEP.OUT.ATTRS = FALSE
)
design$seed <- 500000 + seq_len(nrow(design))

### Simulation
results <- lapply(seq_len(nrow(design)), function(i) {
  d <- design[i, , drop = FALSE]
  set.seed(d$seed)

  terms <- sprintf("log_hazard(t, delta = %s)", d$delta)
  if (d$covariate_setting == "none") {
    fit_formula <- Surv(time, status) ~ offset(offset) | tend
    covar_spec <- NULL
  } else {
    fit_formula <- Surv(time, status) ~ x1 + x2 + offset(offset) | tend
    terms <- c(terms, "0.7 * x1", "0.4 * x2")
    covar_spec <- covariate_spec
  }
  sim_formula <- reformulate(terms)

  out <- tryCatch(
    {
      sim_data <- sim_pam(
        n = d$n,
        formula = sim_formula,
        covariate_spec = covar_spec,
        censoring_rate = d$censoring_rate,
        admin_time = d$admin_time,
        cut = seq(0, d$admin_time, by = d$sim_interval)
      )

      fit_normal <- pemtree_splines(
        formula = fit_formula,
        data = sim_data,
        cut = seq(0, d$admin_time, by = d$ped_interval),
        min_events = 20,
        cluster_id = NULL,
        alpha = 0.05,
        maxdepth = 2
      )

      fit_cluster <- pemtree_splines(
        formula = fit_formula,
        data = sim_data,
        cut = seq(0, d$admin_time, by = d$ped_interval),
        min_events = 20,
        cluster_id = "id",
        alpha = 0.05,
        maxdepth = 2
      )

      fits <- list(normal = fit_normal, cluster = fit_cluster)
      tree_out <- Map(
        function(fit, method) {
          n_terminal <- length(nodeids(fit$tree, terminal = TRUE))
          split_detected <- n_terminal > 1

          tau_hat <- if (split_detected) {
            root <- node_party(fit$tree)
            sp <- split_node(root)
            sp$breaks[1]
          } else {
            NA_real_
          }

          cbind(
            d,
            method = method,
            success = TRUE,
            n_terminal_nodes = n_terminal,
            split_detected = split_detected,
            tau_hat = tau_hat,
            observed_censoring = mean(sim_data$status == 0),
            error = NA_character_
          )
        }, fits, names(fits)
      )
      out <- do.call(rbind, tree_out)
    },
    error = function(e) {
      cbind(
        d,
        method = c("normal", "cluster"),
        success = FALSE,
        n_terminal_nodes = NA_integer_,
        split_detected = NA,
        tau_hat = NA_real_,
        observed_censoring = NA_real_,
        error = conditionMessage(e)
      )
    }
  )
})

results <- do.call(rbind, results)
rownames(results) <- as.character(seq_len(nrow(results)))
head(results)

## Analysis of results
all.equal(sum(results$success), nrow(results))
data.table::setDT(results)

# Alpha
alpha_results <- results[success == TRUE & delta == 0, .(
  alpha_hat = mean(split_detected),
  mcse = sqrt(mean(split_detected) * (1 - mean(split_detected)) / .N)
), by = .(n, covariate_setting, method)]
setorder(alpha_results, n, covariate_setting, method)
alpha_results

# Power
power_results <- results[success == TRUE & delta == 0.5, .(
  power_hat = mean(split_detected),
  mcse = sqrt(mean(split_detected) * (1 - mean(split_detected)) / .N)
), by = .(n, covariate_setting, method)]
setorder(power_results, n, covariate_setting, method)
power_results
