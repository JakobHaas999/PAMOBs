run_one <- function(pars, dgp, fit_model, evaluate) {
  result <- tryCatch(
    {
      set.seed(pars$seed)

      dat <- dgp(pars)
      model <- fit_model(dat, pars)
      metrics <- evaluate(model, dat, pars)

      if (!is.data.frame(metrics) || nrow(metrics) != 1L) {
        stop("evaluate has to give a data.frame with exactly one row.")
      }

      data.frame(success = TRUE, error = NA_character_, metrics, check.names = FALSE)
    },
    error = function(e) {
      data.frame(success = FALSE, error = conditionMessage(e))
    }
  )

  result
}