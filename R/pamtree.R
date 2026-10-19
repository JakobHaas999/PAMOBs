pamtree <- function(
  formula,
  data,
  cut = NULL,
  min_events = 20,
  sp = NULL,
  gam_method = "REML",
  cluster_id = NULL,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)

  # asserts
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_numeric(cut, any.missing = FALSE, null.ok = TRUE, lower = 0)
  checkmate::assert_number(sp, lower = 0, null.ok = TRUE)
  checkmate::assert_string(gam_method)

  # Construct formulas
  f <- Formula::Formula(formula)
  if (length(f)[[2]] != 2) stop("A two part formula with partitioning variables is required")

  node_formula <- formula(f, lhs = 1, rhs = 1)
  node_formula <- update(node_formula, ped_status ~ .)

  data_vars <- intersect(all.vars(f), colnames(data))
  response <- f[[2]]
  ped_terms <- setdiff(data_vars, all.vars(response))
  ped_formula <- reformulate(ped_terms, deparse(response))


  # Transform to PED data
  ped_args <- list(formula = ped_formula, data = data)
  if (!is.null(cut)) ped_args$cut <- cut
  ped <- do.call(pammtools::as_ped, ped_args)
  ped$.row_id <- seq_len(nrow(ped))

  fit_gam <- function(data, sp = NULL) {
    gam_args <- list(formula = node_formula, family = poisson(link = "log"), data = data, method = gam_method)
    if (!is.null(sp)) gam_args$sp <- sp
    do.call(mgcv::gam, gam_args)
  }

  template_fit <- tryCatch(
    suppressWarnings(fit_gam(ped)),
    error = function(e) {
      stop("Could not fit root GAM template:", conditionMessage(e))
    }
  )

  n_sp <- length(template_fit$sp)
  fixed_sp <- NULL
  if (!is.null(sp)) {
    if (!is.numeric(sp) || any(!is.finite(sp)) || any(sp < 0)) {
      stop("sp must be a finite non-negative numeric vectors")
    }

    fixed_sp <- if (length(sp) == 1) rep(sp, n_sp) else sp
    if (length(fixed_sp) != n_sp) {
      stop(
        "sp must have length 1 or length equal to the number of smoothing parameters: ",
        n_sp
      )
    }
  }

  template_fit <- if (is.null(fixed_sp)) {
    template_fit
  } else {
    tryCatch(
      suppressWarnings(fit_gam(ped, sp = fixed_sp)),
      error = function(e) {
        stop("Could not fit root GAM template:", conditionMessage(e))
      }
    )
  }

  coef_template <- coef(template_fit)

  make_fit_fun <- function(ped, node_formula, min_events, coef_template, gam_method, fixed_sp) {
    force(ped)
    force(node_formula)
    force(min_events)
    force(coef_template)
    force(gam_method)
    force(fixed_sp)

    function(y, x, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
      row_ids <- as.integer(y)
      node_data <- ped[row_ids, , drop = FALSE]

      invalid_fit <- function() {
        scores <- NULL

        if (estfun) {
          scores <- matrix(0, nrow = length(nrow(node_data)), ncol = length(coef_template))
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
        family = poisson(link = "log"),
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
          any(!is.finite(coef(fit)))
      ) {
        return(invalid_fit())
      }

      X <- predict(fit, type = "lpmatrix")
      colnames(X) <- names(coef_template)
      raw_scores <- sweep(
        X, 1, node_data$ped_status - fitted(fit),
        FUN = "*"
      )
      scores <- sweep(
        raw_scores, 2, colMeans(raw_scores),
        FUN = "-"
      )

      list(
        coefficients = coef(fit),
        objfun = -as.numeric(logLik(fit)),
        estfun = if (estfun) scores else NULL,
        object = if (object) fit else NULL
      )
    }
  }

  fit_fun <- make_fit_fun(
    ped, node_formula, min_events,
    coef_template, gam_method, fixed_sp
  )

  partition_f <- formula(f, lhs = 0, rhs = 2)
  partition_terms <- attr(terms(partition_f), "term.labels")
  mob_formula <- Formula::Formula(
    as.formula(paste0(".row_id ~ 1 | ", paste(partition_terms, collapse = " + ")), env = environment())
  )

  mob_args <- list(formula = mob_formula, data = ped, fit = fit_fun, control = control)
  if (!is.null(cluster_id)) {
    if (!(cluster_id %in% colnames(ped))) stop("cluster_id not found in PED data.")
    mob_args$cluster <- ped[[cluster_id]]
  }

  tree <- do.call(partykit::mob, mob_args)

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
    class = c("pamtree", "pemtree")
  )
}


## Example call
# Simulated some data with a structural change
n <- 2000
shape <- 1.5
rate <- 0.05

tau <- 5
delta <- log(5)

log_hazard <- ~ log(rate) + log(shape) + (shape - 1) * log(t) + delta * (t >= tau)
data <- data.frame(id = 1:n)

admin_time <- 10
censoring_rate <- 0.3

set.seed(1827)
sim_data <- sim_surv(
  formula = log_hazard,
  data = data,
  breaks = tau,
  admin_time = admin_time,
  censoring_rate = censoring_rate
)

# Fit the tree
formula <- Surv(time, status) ~ s(tend, bs = "ps") + offset(offset) | tend
tree_fit <- pamtree(
  formula = formula,
  data = sim_data,
  cut = seq(0, 10, by = 0.1),
  maxdepth = 2,
  alpha = 0.05
)
