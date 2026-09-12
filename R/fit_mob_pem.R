make_pem_mob_fit <- function(surv_data,
                             pem_formula,
                             cut) {
  force(surv_data)
  force(pem_formula)
  force(cut)

  function(y, x = NULL,
           start = NULL,
           weights = NULL,
           offset = NULL,
           ...,
           estfun = FALSE,
           object = FALSE) {
    ids <- as.integer(y)
    node_data <- surv_data[surv_data$id %in% ids, , drop = FALSE]

    pem <- fit_pem(formula = pem_formula, data = node_data, cut = cut)

    fit <- pem$fit
    ped <- pem$ped

    scores <- NULL

    if (estfun) {
      scores_row <- extract_scores_pem(fit)
      scores <- aggregate_scores_pem(scores_row, ped)
      scores <- scores[match(ids, as.integer(rownames(scores))), , drop = FALSE]
    }


    list(
      coefficients = coef(fit),
      objfun = -as.numeric(logLik(fit)),
      estfun = scores,
      object = if (object) pem else NULL
    )
  }
}

fit_mob_pem <- function(pem_formula,
                        z_formula,
                        data,
                        cut,
                        ...) {
  z_variables <- attr(terms(z_formula), "term.labels")
  data$id <- seq_len(nrow(data))
  mob_data <- data[, c("id", z_variables)]

  mob_formula <- reformulate(
    termlabels = paste0("1 | ", paste(z_variables, collapse = " + ")),
    response = "id"
  )

  pem_fit <- make_pem_mob_fit(
    surv_data = data,
    pem_formula = pem_formula,
    cut = cut
  )

  partykit::mob(
    formula = mob_formula,
    data = mob_data,
    fit = pem_fit,
    ...
  )
}
