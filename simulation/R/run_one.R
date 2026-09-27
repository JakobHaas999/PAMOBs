extract_time_split <- function(fit) {
  root <- partykit::node_party(fit$tree)
  split <- partykit::split_node(root)

  if (is.null(split)) {
    return(NA_real_)
  }

  partykit::breaks_split(split)[[1]]
}

run_one_level_shift <- function(
  n,
  delta,
  tau,
  beta0,
  censoring_rate,
  rep,
  seed,
  ped_interval = 0.1,
  alpha = 0.05
) {
  cat(sprintf("\nRepetition %s\n", rep))
  start_time <- proc.time()[["elapsed"]]

  result <- tryCatch(
      {
        set.seed(seed)
        censoring_argument <- if (is.null(censoring_rate) || censoring_rate == 0) {
          NULL
        } else {
          censoring_rate
        }
        data <- scenario_level_shift(
          n = n,
          delta = delta,
          tau = tau,
          beta0 = beta0,
          censoring_rate = censoring_argument
        )
        max_time <- max(data$time)
        ped_cut <- unique(c(
          seq(0, max_time, by = ped_interval),
          max_time
        ))
        fit <- fit_temporal_pemtree(
          data = data,
          cut = ped_cut,
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
    split_detected <- !is.na(result$tau_hat)

    data.frame(
      scenario = "level_shift",
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
      tau_error = if (split_detected && delta != 0) {
        abs(result$tau_hat - tau)
      } else {
        NA_real_
      },
      n_terminal_nodes = result$n_terminal_nodes,
      observed_censoring = result$observed_censoring,
      runtime_sec = runtime,
      error = result$error,
      alpha = alpha
    )
}
