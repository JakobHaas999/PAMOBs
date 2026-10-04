### setup
# packages
source("setup.R")
# simulation functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
# pemtree functions
source("R/pemtree.R")
source("R/pemtree_splines.R")


##### Functions -------------------------------------------------------------
## Data-Generating Process
level_shift <- function(n, delta, tau, beta0, censoring_rate, cut) {
  log_hazard <- function(t) {
    beta0 + delta * (t >= tau)
  }

  data <- sim_pam(
    n = n,
    formula = ~ log_hazard(t),
    covariate_spec = NULL,
    censoring_rate = censoring_rate,
    admin_time = max(cut),
    cut = cut
  )

  list(sim_data = data, FUN = log_hazard)
}

## Fit-function
tree_fit <- function(method = c("normal", "splines"), data, cut, min_events = 20, ...) {
  method <- match.arg(method)
  treefun <- if (method == "splines") pemtree_splines else pemtree

  treefun(
    formula = survival::Surv(time, status) ~ 1 + offset(offset) | tend,
    data = data,
    cut = cut,
    min_events = min_events,
    ...
  )
}

## Function to extract time split
extract_split <- function(fit) {
  root <- node_party(fit$tree)
  split <- split_node(root)

  if (is.null(split)) {
    return(NA_real_)
  }

  breaks_split(split)[[1]]
}

## Function to run one simulation
run_one_sim <- function(n,
                        delta,
                        tau,
                        beta0,
                        censoring_rate,
                        seed,
                        rep,
                        method = c("normal", "splines"),
                        cut = seq(0, 10, by = 0.1),
                        min_events = 20,
                        maxdepth = 2,
                        alpha = 0.05) {
  method <- match.arg(method)
  runtime <- NA_real_

  result <- tryCatch(
    {
      set.seed(seed)
      data <- level_shift(n, delta, tau, beta0, censoring_rate, cut = cut)$sim_data
      start_time <- proc.time()[["elapsed"]]
      fit <- tree_fit(method,
        data = data,
        cut = cut,
        min_events = min_events,
        maxdepth = maxdepth,
        alpha = alpha
      )
      runtime <- proc.time()[["elapsed"]] - start_time

      tau_hat <- extract_split(fit)
      n_terminal_nodes <- length(nodeids(fit$tree, terminal = TRUE))

      list(
        success = TRUE, tau_hat = tau_hat, n_terminal_nodes = n_terminal_nodes,
        obs_censoring = mean(data$status == 0), error = NA_character_
      )
    },
    error = function(e) {
      list(
        success = FALSE, tau_hat = NA_real_, n_terminal_nodes = NA_integer_,
        obs_censoring = NA_real_, error = conditionMessage(e)
      )
    }
  )

  split_detected <- if (result$success) {
    !is.na(result$tau_hat)
  } else {
    NA
  }

  data.table::data.table(
    method = method,
    n = n,
    delta = delta,
    tau = tau,
    beta0 = beta0,
    seed = seed,
    success = result$success,
    split_detected = split_detected,
    tau_hat = result$tau_hat,
    n_terminal_nodes = result$n_terminal_nodes,
    obs_censoring = result$obs_censoring,
    error = result$error,
    runtime = runtime
  )
}


## Benchmark study
design <- expand.grid(
  n = 1000,
  delta = c(0, 0.5, 1),
  tau = 4,
  beta0 = -2,
  censoring_rate = 0.2,
  method = c("normal", "splines"),
  rep = seq_len(20),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design$seed <- NA_integer_
design$seed[design$method == "normal"] <- 10L + seq_len(sum(design$method == "normal"))
design$seed[design$method == "splines"] <- 10L + seq_len(sum(design$method == "splines"))

results <- data.table::rbindlist(
  lapply(seq_len(nrow(design)), function(i) {
    do.call(run_one_sim, as.list(design[i, ]))
  })
)

results[, .(
  n_runs = .N,
  n_errors = sum(!success),
  split_rate = mean(split_detected, na.rm = TRUE),
  mean_tau_hat = mean(tau_hat, na.rm = TRUE),
  median_runtime = median(runtime, na.rm = TRUE)
), by = c("method", "delta")]
