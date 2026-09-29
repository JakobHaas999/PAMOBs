# Fit a Piecewise Exponential Model Tree
#
# This function combines piecewise exponential models (PEMs) with
# model-based recursive partitioning. The supplied survival data are first
# transformed into piecewise exponential data (PED). A Poisson GLM tree is
# then fitted to the PED representation using `partykit::glmtree()`.
#
# The formula follows the model-based recursive partitioning structure
#
#   Surv(time, status) ~ model terms | partitioning variables
#
# where the terms on the left-hand side of `|` define the node-specific
# PEM and the variables on the right-hand side define potential splitting
# variables.
#
# The fitted tree, PED representation, model formulas, cut points, and
# additional information required for subsequent methods are stored in
# an object of class "pemtree".
pemtree <- function(
  formula,
  data,
  cut,
  min_events = 10,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_numeric(
    cut,
    any.missing = FALSE,
    lower = 0
  )
  checkmate::assert_count(min_events, positive = TRUE)

  # construct PED formula
  f <- Formula::Formula(formula)
  data_vars <- intersect(all.vars(f), colnames(data))
  response <- f[[2]]
  vars_ped <- setdiff(data_vars, all.vars(response))
  ped_formula <- reformulate(
    termlabels = vars_ped,
    response = deparse(response)
  )

  # create PED data
  ped <- pammtools::as_ped(
    formula = ped_formula,
    data = data,
    cut = cut
  )

  # Custom robust fit function for mob()
  make_poisson_model_robust <- function(min_events) {
    force(min_events)
    function(y, x = NULL,
             start = NULL,
             weights = NULL,
             offset = NULL,
             ...,
             estfun = FALSE,
             object = FALSE) {
      # break criteria
      sum_events <- sum(y)
      if (sum_events < min_events) {
        return(list(
          coefficients = setNames(numeric(ncol(x)), colnames(x)),
          objfun = Inf,
          estfun = NULL,
          object = NULL
        ))
      }

      fit <- glm.fit(
        x = x,
        y = y,
        weights = weights,
        start = start,
        offset = offset,
        family = poisson(link = "log")
      )

      # Reference partykit:::glmfit
      objfun <- fit$aic / 2 - fit$rank

      scores <- NULL
      if (estfun) {
        wres <- as.vector(fit$residuals) * fit$weights
        scores <- wres * x[, !is.na(fit$coefficients), drop = FALSE]
      }

      if (object) {
        class(fit) <- c("glm", "lm")
      }

      return(list(
        coefficients = coef(fit),
        objfun = objfun,
        estfun = scores,
        object = if (object) fit else NULL
      ))
    }
  }
  fit_fun <- make_poisson_model_robust(min_events)

  # Construct mob formula
  mob_formula <- update(f, ped_status ~ .)

  tree <- partykit::mob(
    formula = mob_formula,
    data = ped,
    fit = fit_fun,
    control = control
  )

  structure(list(
    tree = tree,
    formula = formula,
    mob_formula = mob_formula,
    ped_formula = ped_formula,
    ped = ped,
    nobs = nrow(data),
    nobs_ped = nrow(ped),
    cut = cut,
    min_events = min_events,
    call = cl
  ), class = "pemtree")
}

# Print a Piecewise Exponential Model Tree
#
# Prints a compact representation of a fitted "pemtree" object, including
# the original function call and the underlying generalized linear model
# tree.
print.pemtree <- function(x, ...) {
  cat("Piecewise Exponential Model Tree\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nFitted tree:\n")
  print(x$tree, ...)

  invisible(x)
}
