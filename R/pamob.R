### Fit Model-Based Trees with Piecewise Exponential Additive Models in the nodes
#
# This pipeline combines Piecewise Exponential Additive Models (PAMs) with
# Model-Based Recursive Partitioning (MOB)
# # Example call
# formula <- ~ log_hazard(t) + 0.7 * x1 + sqrt(x2) + 0.8 * (x3 == 1) * x1
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
# res <- pamob(
#   formula = Surv(time, status) ~ s(tend, bs = "ps") + x1 + s(x2, bs = "ps") + offset(offset) | x3,
#   data = sim_data,
#   cut = cut,
#   min_events = min_events,
#   alpha = 0.05,
#   maxdepth = 2
# )

pamob <- function(
  formula,
  data,
  cut,
  min_events = 20,
  gam_method = "REML",
  sp = NULL,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)
  data <- as.data.frame(data)
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_numeric(cut, lower = 0, any.missing = FALSE)
  # TODO check cut more
  checkmate::assert_string(gam_method, null.ok = TRUE)

  f <- Formula::Formula(formula)
  if (length(f)[[2]] != 2) {
    stop("A two part formula with partitioning variables is required")
  }

  node_formula <- formula(f, lhs = 1, rhs = 1)
  node_formula <- update(node_formula, ped_status ~ .)

  data_vars <- intersect(colnames(data), all.vars(f))
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


  person_data <- data[, c("id", partition_vars), drop = FALSE]
  if (anyDuplicated(person_data$id)) {
    stop("Partitioning variables must be constant within persons for person-wise scores.")
  }

  fit_gam <- function(data, sp = NULL) {
    gam_args <- list(
      formula = node_formula,
      family = stats::poisson(link = "log"),
      data = data,
      method = gam_method
    )

    if (!is.null(sp)) {
      gam_args$sp <- sp
    }

    do.call(mgcv::gam, gam_args)
  }

  template_fit_reml <- tryCatch(
    suppressWarnings(fit_gam(ped)),
    error = function(e) {
      stop("Could not fit root GAM template", conditionMessage(e))
    }
  )

  n_sp <- length(template_fit_reml$sp)

  fixed_sp <- NULL
  if (!is.null(sp)) {
    if (!is.numeric(sp) || any(!is.finite(sp)) || any(sp < 0)) {
      stop("sp must be a finite non-negative numeric vectors")
    }

    fixed_sp <- if (length(sp) == 1L) rep(sp, n_sp) else sp
    if (length(fixed_sp) != n_sp) {
      stop(
        "sp must have length 1 or length equal to the number of smoothing parameters: ",
        n_sp
      )
    }
  }

  template_fit <- if (is.null(fixed_sp)) {
    template_fit_reml
  } else {
    tryCatch(
      suppressWarnings(fit_gam(ped, sp = fixed_sp)),
      error = function(e) {
        stop("Could not fit root GAM template with fixed sp: ", conditionMessage(e))
      }
    )
  }

  coef_template <- coef(template_fit)

  make_fit_fun <- function(ped, node_formula, min_events,
                           coef_template, gam_method, fixed_sp) {
    force(ped)
    force(node_formula)
    force(min_events)
    force(coef_template)
    force(gam_method)
    force(fixed_sp)

    function(y, x, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
      person_ids <- as.character(y)
      node_data <- ped[as.character(ped$id) %in% person_ids, , drop = FALSE]

      invalid_fit <- function() {
        scores <- NULL

        if (estfun) {
          scores <- matrix(0, length(person_ids), ncol = length(coef_template))
          colnames(scores) <- names(coef_template)
        }
        list(
          coefficients = coef_template * 0,
          objfun = Inf,
          estfun = scores,
          object = NULL
        )
      }

      if (!nrow(node_data) || sum(node_data$ped_status) < min_events) {
        return(invalid_fit())
      }

      gam_args <- list(
        formula = node_formula,
        family = stats::poisson(link = "log"),
        data = node_data,
        method = gam_method
      )

      if (!is.null(fixed_sp)) {
        gam_args$sp <- fixed_sp
      }

      fit <- tryCatch(
        suppressWarnings(do.call(mgcv::gam, gam_args)),
        error = function(e) NULL
      )

      if (
        is.null(fit) ||
          !inherits(fit, "gam") ||
          !isTRUE(fit$converged) ||
          any(!is.finite(stats::coef(fit)))
      ) {
        return(invalid_fit())
      }

      X <- predict(fit, type = "lpmatrix")
      colnames(X) <- names(coef(fit))

      raw_scores <- sweep(
        X, 1, node_data$ped_status - fitted(fit),
        FUN = "*"
      )
      scores <- rowsum(raw_scores, group = node_data$id)
      scores <- scores[match(person_ids, rownames(scores)), , drop = FALSE]

      scores[is.na(scores)] <- 0

      scores <- sweep(
        scores, 2, colMeans(scores),
        FUN = "-"
      )
      colnames(scores) <- names(coef(fit))

      list(
        coefficients = coef(fit),
        objfun = -as.numeric(stats::logLik(fit)),
        estfun = if (estfun) scores else NULL,
        object = if (object) fit else NULL
      )
    }
  }

  fit_fun <- make_fit_fun(ped, node_formula, min_events, coef_template, gam_method, fixed_sp)

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
    mob_formula = mob_formula,
    ped_formula = ped_formula,
    ped = ped,
    person_data = person_data,
    nobs = nrow(data),
    nobs_ped = nrow(ped),
    cut = cut,
    min_events = min_events,
    gam_method = gam_method,
    sp = fixed_sp,
    call = cl
  )
  class(out) <- c("pamob", "pemob")
  out
}
