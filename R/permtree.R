permtree <- function(formula, data, n_perm = 999L, alpha = 0.05) {
  f <- Formula::Formula(formula)

  # Build PED formula
  data_vars <- intersect(colnames(data), all.vars(formula))
  response <- f[[2]]
  ped_terms <- setdiff(data_vars, all.vars(response))
  ped_formula <- reformulate(
    termlabels = ped_terms, response = deparse(response)
  )

  # Transform to PED
  ped <- pammtools::as_ped(
    data = data,
    formula = ped_formula,
    id = "id"
  )

  # Build the node-model formula
  model_formula <- formula(f, lhs = 1, rhs = 1)
  model_formula[[2L]] <- quote(ped_status)

  # Fit the Poisson model
  fit <- glm(
    formula = model_formula,
    family = poisson(link = "log"),
    data = ped,
    na.action = na.fail
  )

  # Extract split variable
  split_formula <- formula(f, lhs = 0, rhs = 2)
  split_vars <- all.vars(split_formula)
  if (length(split_vars) != 1) stop("Exactly one split variable is supported.")
  split_var <- split_vars[[1]]
  z <- droplevels(factor(ped[[split_var]]))

  test <- strucchange::sctest(
    fit,
    order.by = z, functional = "LMuo"
  )

  statistic <- unname(test$statistic)

  person_ids <- unique(ped$id)
  person_index <- match(ped$id, person_ids)
  first_row <- match(person_ids, ped$id)

  z_person <- z[first_row]

  permutation_statictics <- replicate(n_perm, {
    z_perm_person <- sample(z_person)
    z_perm_row <- z_perm_person[person_index]

    test <- strucchange::sctest(
      fit,
      order.by = z_perm_row, functional = "LMuo"
    )
    unname(test$statistic)
  })

  p_val <- (1 + sum(permutation_statictics >= statistic)) / (n_perm + 1)

  if (p_val <= alpha) {
    data_left <- data[data[[split_var]] == 0, ]
    data_right <- data[data[[split_var]] == 1, ]

    ped_left <- pammtools::as_ped(
      data = data_left, formula = ped_formula, id = "id"
    )
    ped_right <- pammtools::as_ped(
      data = data_right, formula = ped_formula, id = "id"
    )

    fit_left <- glm(
      formula = model_formula,
      family = poisson(link = "log"),
      data = ped_left
    )
    fit_right <- glm(
      formula = model_formula,
      family = poisson(link = "log"),
      data = ped_right
    )

    return(list(
      models = list(fit_left, fit_right),
      terminal_nodes = 2,
      p_val = p_val
    ))
  } else {
    return(list(
      models = list(fit),
      terminal_nodes = 1,
      p_val = p_val
    ))
  }
}

set.seed(123)
rate <- 0.05
shape <- 1.5
formula <- ~ log(rate) + log(shape) + (shape - 1) * log(t) + 0.7 * x + 0.6 * x * (z == 1)

n_obs <- 1000
n_pilot <- 10000

censoring_rate <- 0.3
admin_time <- 10

event_times <- sim_event_times(
  formula = formula,
  data = data.frame(x = rnorm(n_pilot), z = rbinom(n_pilot, 1, 0.5), id = seq_len(n_pilot))
)
lambda_c <- calibrate_censoring(event_times,
  censoring_rate = censoring_rate,
  admin_time = admin_time
)

sim_data <- sim_surv(
  formula = formula,
  data = data.frame(x = rnorm(n_obs), z = rbinom(n_obs, 1, 0.5), id = seq_len(n_obs)),
  lambda_c = lambda_c,
  admin_time = admin_time
)

res <- permtree(
  formula = Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | z,
  data = sim_data,
  n_perm = 100L
)
