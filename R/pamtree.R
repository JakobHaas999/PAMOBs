pamtree <- function(
  formula,
  data,
  cut = NULL,
  min_events,
  gam_method = "REML",
  sp,
  cluster_id = NULL,
  ...
) {
  cl <- match.call()
  control <- partykit::mob_control(...)

  # Asserts
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_numeric(cut, lower = 0, any.missing = FALSE)
  checkmate::assert_count(min_events, positive = TRUE)
  checkmate::assert_string(cluster_id, null.ok = TRUE)

  # Construct PED formula
  f <- Formula::Formula(formula)
  if (length(f)[[2]] != 2) {
    stop("A two part formula is partitioning variables is required")
  }

  node_formula <- formula(f, lhs = 1, rhs = 1)
  node_formula <- update(node_formula, ped_status ~ .)


  data_vars <- intersect(all.vars(f), colnames(data))
  response <- f[[2]]
  ped_vars <- setdiff(data_vars, all.vars(response))
  ped_formula <- reformulate(
    termlabels = ped_vars,
    response = deparse(response)
  )

  ped_args <- list(formula = ped_formula, data = data)
  if (!is.null(cut)) ped_args$cut <- cut
  ped <- do.call(pammtools::as_ped, ped_args)
  ped$.row_id <- seq_len(nrow(ped))

  make_node_pam_fit <- function(ped, node_formula, min_events, gam_method, sp) {
    force(ped)
    force(node_formula)
    force(min_events)
    force(gam_method)
    force(sp)

    function(y, x = NULL,
             start = NULL,
             weigths = NULL,
             offset = NULL,
             ...,
             estfun = FALSE,
             object = FALSE) {
      row_ids <- as.integer(y)
      node_data <- ped[row_ids, , drop = FALSE]

      if (sum(node_data$ped_status) < min_events) {}
    }
  }
}
