# NOTE:
# PED rows from the same subject are not independent.
# This analysis is currently exploratory and is used to assess whether
# a temporal structural break is visible to MOB at all.

#' Fit a temporal PEM-MOB
#'
#' Transforms survival data into piecewise exponential data (PED) and fits
#' a model-based recursive partitioning tree using time (`tend`) as the
#' partitioning variable. Within each node, a Poisson model corresponding
#' to a piecewise exponential model is fitted.
#'
#' @param formula Survival formula of the form Surv(time, status) ~ covariates.
#' @param data Data frame containing survival outcome and covariates.
#' @param cut Cut points used for the PED transformation.
#' @param ... Additional arguments passed to partykit::glmtree().
#'
#' @return An object of class "temporalPemob".
temporal_pemob <- function(formula, data, cut, ...) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_numeric(cut, any.missing = FALSE)

  covariates <- attr(terms(formula), "term.labels")

  # Transform survival data into piecewise exponential data
  ped <- pammtools::as_ped(
    formula = formula,
    data = data,
    cut = cut
  )

  # Use the original covariates in the node model and time as the
  # partitioning variable. Without covariates, fit an intercept-only model.
  tree_formula <- if (length(covariates)) {
    as.formula(paste0(
      "ped_status ~ ",
      paste(covariates, collapse = " + "),
      " + offset(offset) | tend"
    ))
  } else {
    ped_status ~ 1 + offset(offset) | tend
  }

  fit <- partykit::glmtree(
    formula = tree_formula,
    family = poisson(link = "log"),
    data = ped,
    ...
  )

  out <- structure(list(
    tree = fit,
    formula = formula,
    tree_formula = tree_formula,
    data = data,
    ped = ped,
    cut = cut
  ), class = "temporalPemob")

  out
}

#' Print a temporal PEM-MOB
#'
#' Prints the underlying model-based recursive partitioning tree.
print.temporalPemob <- function(x, ...) {
  print(x$tree, ...)
}


#' Predict from a temporal PEM-MOB

#'

#' Assigns each row of `newdata` to a terminal node of the fitted tree.
#' Predictions can either return the corresponding node ID or the estimated
#' hazard. Hazard predictions are obtained from the node-specific Poisson
#' model with an offset of zero, corresponding to unit exposure time.
#'
#' @param object A fitted "temporalPemob" object.
#' @param newdata Data frame containing `tend` and all covariates required
#'   by the fitted node model.
#' @param type Either "node" or "hazard".
#' @param ... Additional arguments.
#'
#' @return A vector containing node IDs or estimated hazards.
predict.temporalPemob <- function(object, newdata, type = c("node", "hazard"), ...) {
  type <- match.arg(type)
  checkmate::assert_data_frame(newdata)
  checkmate::assert_names(colnames(newdata), must.include = "tend")

  nodes <- predict(
    object$tree,
    newdata = newdata,
    type = "node"
  )
  if (type == "node") {
    return(nodes)
  }

  terminal_nodes <- partykit::nodeids(object$tree, terminal = TRUE)
  node_models <- partykit::nodeapply(
    object$tree,
    ids = terminal_nodes,
    FUN = function(node) node$info$object
  )
  names(node_models) <- terminal_nodes

  # Hazard predictions
  hazard <- numeric(nrow(newdata))

  for (node_id in terminal_nodes) {
    idx <- nodes == node_id
    if (!any(idx)) next
    mod <- node_models[[as.character(node_id)]]
    node_data <- newdata[idx, , drop = FALSE]
    node_data$offset <- 0
    hazard[idx] <- predict(
      mod,
      newdata = node_data,
      type = "response",
      ...
    )
  }

  unname(hazard)
}
