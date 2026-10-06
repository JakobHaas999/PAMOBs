construct_tree_formula <- function(covariate_spec, df) {
  model_terms <- c(
    sprintf("bs(tend, degree = 3, df = %s)", df),
    names(covariate_spec),
    "offset(offset)"
  )

  as.formula(
    paste(
      "Surv(time, status) ~",
      paste(model_terms, collapse = " + "),
      "| tend"
    ),
    env = parent.frame()
  )
}

extract_time_split <- function(fit) {
  root <- node_party(fit$tree)
  split <- split_node(root)
  if (is.null(split)) {
    return(NA_real_)
  }
  breaks_split(split)[[1]]
}

run_one_level_shift_splines <- function(
  n,
  delta,
  tau,
  basehaz_fun,
  covariate_setting = c("none", "included"),
  covariate_spec = NULL,
  covariate_effects = NULL,
  df,
  censoring_rate,
  rep,
  seed,
  admin_time,
  ped_interval = 0.1,
  min_events = 25,
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
      censoring_argument <- if (is.null(censoring_rate) || censoring_rate == 0) {
        NULL
      } else {
        censoring_rate
      }

      if (covariate_setting == "none") {
        if (!is.null(covariate_spec) || !is.null(covariate_effects)) {
          stop("NO covariates may be supplied for covariate_setting = 'none'")
        }
      } else {
        if (is.null(covariate_spec) || is.null(covariate_effects)) {
          stop("covariate_spec and covariate_effects must be supplied.")
        }
      }


      data <- scenario_level_shift_splines(
        n = n,
        delta = delta,
        tau = tau,
        basehaz_fun = basehaz_fun,
        censoring_rate = censoring_argument,
        covariate_spec = covariate_spec,
        covariate_effects = covariate_effects,
        admin_time = admin_time,
        sim_interval = ped_interval
      )

      tree_formula <- construct_tree_formula(
        covariate_spec = covariate_spec, df = df
      )

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
      n_terminal_nodes <- length(nodeids(fit$tree, terminal = TRUE))

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
  split_detected <- if (result$success) !is.na(result$tau_hat) else NA

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
    df = df,
    censoring_rate = if (is.null(censoring_rate)) 0 else censoring_rate,
    rep = rep,
    seed = seed,
    success = result$success,
    split_detected = split_detected,
    tau_hat = result$tau_hat,
    observed_censoring = result$observed_censoring,
    runtime_sec = runtime,
    error = result$error,
    admin_time = admin_time,
    ped_interval = ped_interval,
    min_events = min_events,
    alpha = alpha
  )
}
