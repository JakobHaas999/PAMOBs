run_covariate_instability <- function(
  log_hazard,
  design,
  sim_formula,
  fit_formula,
  fit_fun,
  covariate_spec,
  ...
) {
  results <- lapply(seq_len(nrow(design)), function(i) {
    d <- design[i, , drop = FALSE]
    set.seed(d$seed)
    cut <- seq(0, d$admin_time, by = d$sim_interval)
    cut_fit <- seq(0, d$admin_time, by = d$fit_interval)
    sim_formula_i <- if (is.function(sim_formula)) {
      sim_formula(d)
    } else {
      sim_formula
    }

    out <- tryCatch(
      {
        sim_data <- sim_pam(
          n = d$n,
          formula = sim_formula_i,
          covariate_spec = covariate_spec,
          censoring_rate = d$censoring_rate,
          admin_time = d$admin_time,
          cut = cut
        )

        fit <- fit_fun(
          formula = fit_formula,
          data = sim_data,
          cut = cut_fit,
          ...
        )

        n_terminal <- length(partykit::nodeids(fit$tree, terminal = TRUE))
        if (inherits(fit, "pamob")) {
          root_info <- partykit::info_node(partykit::node_party(fit$tree))
          test_failed <- is.null(root_info$test) || anyNA(root_info$test["p.value", ])
          if (test_failed) stop("MOB fluctuation test failed")
        }
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

  do.call(rbind, results)
}
