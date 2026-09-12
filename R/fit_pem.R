### Function for fitting a PEM
#
## Example call
data <- pammtools::tumor[1:100, ] %>%
  select(days, status, charlson_score, sex, age)

fit_pem <- function(formula, data, ...) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)

  ped <- pammtools::as_ped(
    formula,
    data = data,
    ...
  )

  terms_rhs <- attr(terms(formula), "term.labels")

  formula_glm <- reformulate(
    termlabels = c("interval", terms_rhs),
    response = "ped_status"
  )

  fit <- glm(
    formula = formula_glm,
    family = poisson(link = "log"),
    data = ped,
    offset = ped$offset
  )

  list(
    fit = fit,
    ped = ped
  )
}
