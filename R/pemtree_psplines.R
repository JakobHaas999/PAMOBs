pemtree_psplines <- function(
  formula,
  data,
  cut,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)

  f <- Formula::Formula(formula)
  if (length(f)[[2]] != 2) {
    stop("A two part formula with partitioning variables is required")
  }


  node_formula <- formula(f, lhs = 1, rhs = 1)
  node_formula <- update(node_formula, ped_status ~ .)

  data_vars <- intersect(all.vars(f), colnames(data))
  response <- f[[2]]
  ped_formula <- reformulate(
    setdiff(data_vars, all.vars(response)),
    response = deparse(response)
  )

  ped <- pammtools::as_ped(
    formula = ped_formula,
    data = data,
    cut = cut
  )

  ped$.row_id <- seq_len(nrow(ped))

  make_node_pam_fit <- function(ped, node_formula) {
    force(ped)
    force(node_formula)

    function(y, x = NULL,
             start = NULL,
             weights = NULL,
             offset = NULL,
             ...,
             estfun = FALSE,
             object = FALSE) {
      row_ids <- as.integer(y)
      node_data <- ped[row_ids, , drop = FALSE]

      invalid_fit <- function(design) {
        list(
          coefficients = setNames(numeric(ncol(design)), colnames(design)),
          objfun = Inf,
          estfun = NULL,
          object = NULL
        )
      }

      node_data$.mob_weights <- if (is.null(weights)) rep(1, nrow(node_data)) else weights

      fit <- suppressWarnings(mgcv::gam(
        formula = node_formula,
        family = stats::poisson(link = "log"),
        data = node_data,
        weights = .mob_weights,
        method = "REML"
      ))

      x <- predict(fit, type = "lpmatrix")
      colnames(x) <- names(coef(fit))
      if (!isTRUE(fit$converged) ||
        any(!is.finite(coef(fit)))) {
        return(invalid_fit(x))
      }

      w <- node_data$.mob_weights
      wres <- as.numeric(node_data$ped_status - fitted(fit)) * w

      raw_scores <- sweep(x, 1, wres, "*")
      scores <- raw_scores - tcrossprod(
        w / sum(w), colSums(raw_scores)
      )

      list(
        coefficients = coef(fit),
        objfun = -as.numeric(logLik(fit)),
        estfun = if (estfun) scores else NULL,
        object = if (object) fit else NULL
      )
    }
  }

  fit_fun <- make_node_pam_fit(ped, node_formula)

  partition_formula <- formula(f, lhs = 0, rhs = 2)
  partitions_terms <- attr(terms(partition_formula), "term.labels")
  mob_formula <- Formula::Formula(
    as.formula(paste(".row_id ~ 1 |", paste(partitions_terms, collapse = " + ")),
      env = environment(formula)
    )
  )
  tree <- partykit::mob(
    formula = mob_formula,
    data = ped,
    fit = fit_fun
  )

  tree
}


log_hazard <- function(t) {
  ifelse(t < 5,
    -2.4 + 1.4 * t + -0.08 * t^2,
    -2
  )
}

covariate_spec <- list(
  x1 = list(distfun = rnorm, mean = 0, sd = 1),
  x2 = list(distfun = runif, min = -2, max = 2)
)

sim_formula <- ~ log_hazard(t) + 0.4 * x1 + 0.3 * x2
censoring_rate <- 0.2
admin_time <- 10

set.seed(123)
data <- sim_pam(
  n = 1000,
  formula = sim_formula,
  covariate_spec = covariate_spec,
  censoring_rate = censoring_rate,
  admin_time = admin_time,
  cut = seq(0, admin_time, by = 0.05)
)

formula <- Surv(time, status) ~ s(tend, bs = "ps") + x1 + x2 +
  offset(offset) | tend
cut <- seq(0, admin_time, by = 0.1)
min_events <- 20

pemtree_psplines(
  formula, data, cut
)
