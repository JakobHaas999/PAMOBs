pemtree_splines <- function(
  formula,
  data,
  cut,
  min_events = 20,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data, any.missing = FALSE)
  checkmate::assert_numeric(
    cut,
    any.missing = FALSE,
    lower = 0
  )
  checkmate::assert_count(min_events, positive = TRUE)

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

  make_node_pem_fit <- function(ped, node_formula, min_events) {
    force(ped)
    force(node_formula)
    force(min_events)

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

      if (sum(node_data$ped_status) < min_events) {
        design <- suppressWarnings(model.matrix(node_formula, data = node_data))
        return(invalid_fit(design))
      }

      node_data$.mob_weights <- if (is.null(weights)) 1 else weights

      fit <- suppressWarnings(glm(
        node_formula,
        data = node_data,
        family = poisson(link = "log"),
        weights = .mob_weights,
        x = TRUE
      ))

      design <- fit$x
      if (!isTRUE(fit$converged) ||
        fit$rank < ncol(design) ||
        any(!is.finite(coef(fit)))) {
        return(invalid_fit(design))
      }

      list(
        coefficients = coef(fit),
        objfun = -as.numeric(logLik(fit)),
        estfun = if (estfun) sandwich::estfun(fit) else NULL,
        object = if (object) fit else NULL
      )
    }
  }

  fit_fun <- make_node_pem_fit(ped, node_formula, min_events)

  partition_formula <- formula(f, lhs = 0, rhs = 2)
  partition_terms <- attr(terms(partition_formula), "term.labels")
  mob_formula <- Formula::Formula(
    as.formula(paste(".row_id ~ 1 |", paste(partition_terms, collapse = " + ")), env = environment(formula))
  )

  tree <- partykit::mob(
    formula = mob_formula,
    data = ped,
    fit = fit_fun,
    control = control
  )

  structure(
    list(
      tree = tree,
      formula = formula,
      node_formula = node_formula,
      mob_formula = mob_formula,
      ped_formula = ped_formula,
      ped = ped,
      nobs = nrow(data),
      nobs_ped = nrow(ped),
      cut = cut,
      min_events = min_events,
      call = cl
    ),
    class = c("pemtree")
  )
}
