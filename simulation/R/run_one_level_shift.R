run_one_level_shift <- function(
  n,
  delta,
  tau,
  beta0,
  covariate_setting = c("none", "included"),
  covariate_spec = NULL,
  covariate_effects = NULL,
  censoring_rate,
  rep,
  seed,
  admin_time = 10,
  sim_interval = 0.05,
  ped_interval = 0.1,
  min_events = 10,
  alpha = 0.05
) {
  covariate_setting <- match.arg(covariate_setting)
  cat(sprintf("\nRepetition %s\n", rep))
  cat("Parameters:\n")
  cat("n = ", n, "delta = ", delta, "censoring rate = ", censoring_rate, "seed = ", seed)

  start_time <- proc.time()[["elapsed"]]

  result <- tryCatch(
    {
      set.seed(seed)

      if (covariate_setting == "none") {
        if (!is.null(covariate_spec) || !is.null(covariate_effects)) {
          stop("No covariates may be supplied for covariate setting = 'none'")
        }
      } else {
        if (is.null(covariate_spec) || is.null(covariate_effects)) {
          stop("covariate_spec and covariate_effects must be supplied.")
        }
      }

      data <- scenario_level_shift(
        n = n,
        delta = delta,
        tau = tau,
        beta0 = beta0,
        censoring_rate = censoring_rate,
        covariate_spec = covariate_spec,
        covariate_effects = covariate_effects,
        admin_time = admin_time,
        sim_interval = sim_interval
      )

      tree_formula <- construct_tree_formula(covariate_spec)

      fit_cut <- sort(unique(c(
        seq(0, max(data$time), by = ped_interval),
        max(data$time)
      )))

      fit <- pemtree_splines(
        formula = tree_formula,
        data = data,
        cut = fit_cut,
        min_events = min_events,
        maxdepth = 2,
        alpha = alpha
      )

      tau_hat <- extract_time_split(fit)
      n_terminal_nodes <- length(
        partykit::nodeids(fit$tree, terminal = TRUE)
      )

      list(
        success = TRUE,
        tau_hat = tau_hat,
        n_terminal_nodes = n_terminal_nodes,
        observed_censoring = mean(data$status == 0),
        error = NA_character_
      )
    },
    error = function(e) {
      list(
        success = FALSE,
        tau_hat = NA_real_,
        n_terminal_nodes = NA_integer_,
        observed_censoring = NA_real_,
        error = conditionMessage(e)
      )
    }
  )

  runtime <- proc.time()[["elapsed"]] - start_time
  split_detected <- if (result$success) {
    !is.na(result$tau_hat)
  } else {
    NA
  }

  scenario <- if (covariate_setting == "none") {
    "level shift without covariates"
  } else {
    "level shift with covariates"
  }

  data.frame(
    scenario = scenario,
    n = n,
    delta = delta,
    tau = tau,
    beta0 = beta0,
    censoring_rate = if (is.null(censoring_rate)) 0 else censoring_rate,
    rep = rep,
    seed = seed,
    success = result$success,
    split_detected = split_detected,
    tau_hat = result$tau_hat,
    n_terminal_nodes = result$n_terminal_nodes,
    observed_censoring = result$observed_censoring,
    runtime_sec = runtime,
    error = result$error,
    admin_time = admin_time,
    sim_interval = sim_interval,
    ped_interval = ped_interval,
    min_events = min_events,
    alpha = alpha
  )
}
