## TODO add more than depth 2 and add more than binary categories
#' Permutation-Based Tree for Piecewise Exponential Models
#'
#' Fits a Poisson model to piecewise exponential data (PED) and tests parameter
#' instability across categorical, subject-level partitioning variables.
#' Permutation p-values are obtained by jointly permuting partitioning-variable
#' vectors between subjects while preserving each subject's PED rows.
#' After multiple-testing adjustment, the variable with the smallest p-value
#' is selected. If significant, separate models are fitted for its categories.
#' Otherwise, the common model is retained. No further splits are performed.
#'
#' @param formula A two-part formula of the form
#'   Surv(time, status) ~ model_terms + offset(offset) | split_variables.
#' @param data A data frame with one row per subject, containing survival time,
#'   event status, covariates, partitioning variables, and a unique `id` column.
#' @param n_perm Number of random subject-level permutations. Defaults to 999.
#' @param alpha Significance level applied to adjusted p-values. Defaults to 0.05.
#' @param maxdepth Currently only checked to reject values greater than 2;
#'   it does not otherwise control tree growth.
#' @param adjusting_method Multiple-testing adjustment method passed to
#'   stats::p.adjust(). Defaults to "bonferroni".
#'
#' @return A list containing:
#' \describe{
#'   \item{models}{Terminal-node Poisson models, named by category if a split
#'     is performed, or a single model named "root" otherwise.}
#'   \item{root_model}{The common Poisson model fitted before testing.}
#'   \item{terminal_nodes}{Number of terminal-node models.}
#'   \item{p_vals}{Named vector of unadjusted permutation p-values.}
#'   \item{p_adjusted}{Named vector of adjusted p-values.}
#'   \item{split}{Logical indicating whether a split was performed.}
#'   \item{split_var}{Selected partitioning variable, or NULL if no split
#'     was performed.}
#' }
PErMtree <- function(formula, data, n_perm = 999L, alpha = 0.05, cut = NULL,
                     maxdepth = 2, adjusting_method = "bonferroni") {
  if (maxdepth > 2) stop("Not supported yet.")

  f <- Formula::Formula(formula)

  # Build PED formula
  data_vars <- intersect(colnames(data), all.vars(formula))
  response <- f[[2]]
  ped_terms <- setdiff(data_vars, all.vars(response))
  ped_formula <- reformulate(
    termlabels = ped_terms, response = deparse(response)
  )

  # Transform to PED
  ped_args <- list(data = data, formula = ped_formula, id = "id")
  if (!is.null(cut)) ped_args$cut <- cut
  ped <- do.call(pammtools::as_ped, ped_args)

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

  # Extract split variables
  split_formula <- formula(f, lhs = 0, rhs = 2)
  split_vars <- all.vars(split_formula)

  person_ids <- unique(ped$id)
  person_index <- match(ped$id, person_ids)
  first_rows <- match(person_ids, ped$id)

  n_persons <- length(person_ids)
  perm_ind <- replicate(n_perm, sample.int(n_persons))

  # Test each split variable
  p_vals <- setNames(numeric(length(split_vars)), split_vars)
  statistics <- setNames(numeric(length(split_vars)), split_vars)

  for (var in split_vars) {
    z <- factor(ped[[var]])
    z_person <- z[first_rows]

    test <- strucchange::sctest(fit, order.by = z, functional = "LMuo")
    stat <- unname(test$statistic)
    statistics[var] <- stat

    perm_stats <- vapply(seq_len(n_perm), function(b) {
      z_perm_person <- z_person[perm_ind[, b, drop = TRUE]]
      z_perm_row <- z_perm_person[person_index]

      perm_test <- strucchange::sctest(fit, order.by = z_perm_row, functional = "LMuo")
      unname(perm_test$statistic)
    }, numeric(1))

    p_vals[var] <- (1 + sum(perm_stats >= stat)) / (n_perm + 1)
  }

  p_adjusted <- p.adjust(p_vals, method = adjusting_method)
  best_var <- names(p_vals)[which.min(p_vals)]
  split_detected <- unname(p_adjusted[best_var] <= alpha)

  if (split_detected) {
    split_values <- factor(ped[[best_var]])
    groups <- levels(split_values)

    models <- lapply(groups, function(group) {
      node_data <- ped[split_values == group, , drop = FALSE]
      glm(
        formula = model_formula,
        family = poisson(link = "log"),
        data = node_data,
        na.action = na.fail
      )
    })
    names(models) <- groups
  } else {
    models <- list(root = fit)
  }

  list(
    models = models,
    root_model = fit,
    terminal_nodes = length(models),
    p_vals = p_vals,
    p_adjusted = p_adjusted,
    split = split_detected,
    split_var = if (split_detected) best_var else NULL
  )
}

## Example call
# set.seed(123)
# rate <- 0.05
# shape <- 1.5
# formula <- ~ log(rate) + log(shape) + (shape - 1) * log(t) + 0.7 * x + 0.6 * x * (z1 == 1) + 0 * x * (z2 == 1)
#
# n_obs <- 1000
# n_pilot <- 10000
#
# censoring_rate <- 0.3
# admin_time <- 10
#
# event_times <- sim_event_times(
#   formula = formula,
#   data = data.frame(
#     x = rnorm(n_pilot),
#     z1 = rbinom(n_pilot, 1, 0.5),
#     z2 = rbinom(n_pilot, 1, 0.8),
#     id = seq_len(n_pilot)
#   )
# )
# lambda_c <- calibrate_censoring(event_times,
#   censoring_rate = censoring_rate,
#   admin_time = admin_time
# )
#
# sim_data <- sim_surv(
#   formula = formula,
#   data = data.frame(x = rnorm(n_obs), z1 = rbinom(n_obs, 1, 0.5), z2 = rbinom(n_obs, 1, 0.8), id = seq_len(n_obs)),
#   lambda_c = lambda_c,
#   admin_time = admin_time
# )
#
# res <- PErMtree(
#   formula = Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | z1 + z2,
#   data = sim_data,
#   n_perm = 100L
# )
#
# res$p_vals
# lapply(res$models, coef)
