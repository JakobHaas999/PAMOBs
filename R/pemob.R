### Fit Model-Based Trees with Piecewise Exponential Models in the nodes
#
# This pipeline combines Piecewise Exponential Models (PEMs) with
# Model-Based Recursive Partitioning (MOB)
#
## Example call
# formula <- ~ log_hazard(t) + 0.7 * x1 + sqrt(x2) + 0.5 * (x3 == 1) * x1
#
# censoring_rate <- 0.3
# admin_time <- 10
# cut <- seq(0, admin_time, by = 0.01)
#
# set.seed(123)
# sim_data <- sim_pam(
#   n = 1000,
#   formula = formula,
#   covariate_spec = covariate_spec,
#   censoring_rate = censoring_rate,
#   admin_time = admin_time,
#   cut = cut
# )
#
# res <- pemob(
#   formula = Surv(time, status) ~ bs(tend, df = 5) + x1 + bs(x2, df = 5) + offset(offset) | x3,
#   data = sim_data,
#   cut = cut,
#   min_events = min_events,
#   alpha = 0.05,
#   maxdepth = 2
# )
pemob <- function(
  formula,
  data,
  cut,
  min_events = 20,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)
  # TODO asserts
  ###############

  f <- Formula::Formula(formula)
  if (length(f)[[2]] != 2) {
    stop("A two part formula with partitioning variables is required")
  }


  node_formula <- formula(f, lhs = 1, rhs = 1)
  node_formula <- update(node_formula, ped_status ~ .)

  data_vars <- intersect(all.vars(f), colnames(data))
  response <- f[[2]]
  ped_formula <- reformulate(
    termlabels = setdiff(data_vars, all.vars(response)),
    response = deparse(response)
  )

  ped <- pammtools::as_ped(
    formula = ped_formula,
    data = data,
    cut = cut
  )

  partition_formula <- formula(f, lhs = 0, rhs = 2)
  partition_terms <- attr(terms(partition_formula), "term.labels")
  partition_vars <- all.vars(partition_formula)

  person_data <- unique(data[, c("id", partition_vars), drop = FALSE])

  if (anyDuplicated(person_data$id)) {
    stop("Partitioning variables must be constant within persons for person-wise scores.")
  }


  make_fit_fun <- function(ped, node_formula, min_events) {
    force(ped)
    force(node_formula)
    force(min_events)
    function(y, x, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
      person_ids <- y
      node_data <- ped[ped$id %in% person_ids, , drop = FALSE]

      if (sum(node_data$ped_status) < min_events) {
        X <- suppressWarnings(model.matrix(node_formula, data = node_data))
        return(list(
          coefficients = setNames(numeric(ncol(X)), colnames(X)),
          objfun = Inf,
          estfun = NULL,
          object = NULL
        ))
      }


      fit <- glm(
        formula = node_formula,
        family = poisson(link = "log"),
        data = node_data,
        x = TRUE
      )

      X <- fit$x
      if (!isTRUE(fit$converged) ||
        fit$rank < ncol(X) ||
        any(!is.finite(coef(fit)))) {
        return(list(
          coefficients = setNames(numeric(ncol(X)), colnames(X)),
          objfun = Inf,
          estfun = NULL,
          object = NULL
        ))
      }

      scores <- NULL
      if (estfun) {
        scores_row <- sandwich::estfun(fit)
        scores <- rowsum(scores_row, group = node_data$id)
        scores <- scores[match(as.character(person_ids), rownames(scores)), , drop = FALSE]
      }

      list(
        coefficients = coef(fit),
        objfun = -as.numeric(logLik(fit)),
        estfun = scores,
        object = if (object) fit else NULL
      )
    }
  }

  fit_fun <- make_fit_fun(ped, node_formula, min_events)


  mob_formula <- Formula::Formula(
    as.formula(
      paste("id ~ 1 |", paste(partition_terms, collapse = " + ")),
      env = environment(formula)
    )
  )

  tree <- partykit::mob(
    formula = mob_formula,
    data = person_data,
    fit = fit_fun,
    control = control
  )

  out <- list(
    tree = tree,
    formula = formula,
    node_formula = node_formula,
    ped_formula = ped_formula,
    ped = ped,
    nobs = nrow(data),
    nobs_ped = nrow(ped),
    cut = cut,
    min_events = min_events,
    call = cl
  )

  class(out) <- "pemob"
  out
}

## Print method for a pemob object
print.pemob <- function(x, ...) {
  cat("Piecewise Exponential Model Tree\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nFitted tree:\n")
  print(x$tree, ...)

  invisible(x)
}
