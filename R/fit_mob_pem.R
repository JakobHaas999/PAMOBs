### Fit Model-Based Trees with Piecewise Exponential Models in the nodes
#
# This pipeline combines Piecewise Exponential Models (PEMs) with
# Model-Based Recursive Partitioning (MOB).

### Fit a Model-Based Tree with PEMs as node models
#
# The function separates two types of variables:
#
#   pem_formula:
#     Variables that are part of the PEM fitted within each node.
#
#   z_formula:
#     Partitioning variables along which MOB searches for parameter
#     instability.
#
# Example:
#
#   pem_formula = Surv(time, status) ~ charlson_score
#   z_formula   = ~ sex + age
#
# fits the PEM
#
#   log(lambda_ij) = alpha_j + beta * charlson_score_i
#
# and tests whether its parameters are stable with respect to
# sex and age.
#
#
## Example function call
# data <- pammtools::tumor[1:100, ] |>
#   dplyr::select(days, status, charlson_score, sex, age)
#
# tree <- pemob(
#   pem_formula = survival::Surv(days, status) ~ charlson_score,
#   z_formula = ~ sex + age,
#   data = data,
#   cut = c(0, 1000, 2000, max(data$days))
# )
#
# plot(tree)
pemob <- function(pem_formula, z_formula, data, cut, ...) {
  # Arguments:
  # pem_formula: survival formula defining the PEM in each node
  # z_formula: formula defining the MOB partitioning variables
  # data: survival data with one row per individual
  # cut: numeric vector defining the PEM intervals
  # ...: additional arguments passed to partykit::mob()
  # TODO add asserts

  z_vars <- attr(terms(z_formula), "term.labels")

  data$id <- seq_len(nrow(data))

  # Create PED
  pem_data <- prepare_pem(
    formula = pem_formula,
    data = data,
    cut = cut
  )

  ped <- pem_data$ped
  glm_formula <- pem_data$glm_formula

  mob_data <- data[, c("id", z_vars), drop = FALSE]
  mob_formula <- reformulate(
    termlabels = paste0("1 | ", paste(z_vars, collapse = " + ")),
    response = "id"
  )

  pem_fit <- make_pem_mob_fit(ped, glm_formula)
  partykit::mob(
    formula = mob_formula,
    data = mob_data,
    fit = pem_fit
  )
}

#--- Helpers -------------------------------------------
#-------------------------------------------------------

### Prepare the Piecewise Exponential Data (PED) representation
prepare_pem <- function(formula, data, cut) {
  ped <- pammtools::as_ped(
    formula = formula,
    data = data,
    cut = cut
  )

  terms_rhs <- attr(terms(formula), "term.labels")
  glm_formula <- reformulate(
    termlabels = c("interval", terms_rhs, "offset(offset)"),
    response = "ped_status"
  )

  list(
    ped = ped,
    glm_formula = glm_formula
  )
}

### Create the node fitting function required by partykit::mob
#
# This function returns a fitting function that MOB can call repeatedly.
# The full PED data set is stored in the closure and is therefore not
# reconstructed within each node.
#
# For every MOB node:
#   1. Determine which individuals belong to the node.
#   2. Select all PED rows belonging to these individuals.
#   3. Fit the node-specific Poisson PEM.
#   4. Calculate score contributions if requested by MOB.
#   5. Aggregate interval-level scores to subject-level scores.
make_pem_mob_fit <- function(ped, glm_formula) {
  force(ped)
  force(glm_formula)

  function(y, x = NULL,
           start = NULL,
           weights = NULL,
           offset = NULL,
           ...,
           estfun = FALSE,
           object = FALSE) {
    ids <- as.integer(y)
    node_data <- ped[ped$id %in% ids, , drop = FALSE]

    # Fit the Poisson model
    fit <- glm(
      formula = glm_formula,
      family = poisson(link = "log"),
      data = node_data
    )

    scores <- NULL

    if (estfun) {
      scores_row <- sandwich::estfun(fit)
      scores <- rowsum(scores_row, group = node_data$id)
      scores <- scores[match(ids, as.integer(rownames(scores))), , drop = FALSE]
    }

    list(
      coefficients = coef(fit),
      objfun = -as.numeric(logLik(fit)),
      estfun = scores,
      object = if (object) fit else NULL
    )
  }
}
