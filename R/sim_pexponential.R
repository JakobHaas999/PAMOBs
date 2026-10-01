sim_pexponential <- function(formula, data, cut) {
  checkmate::assert_formula(formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_numeric(
    cut,
    any.missing = FALSE,
    lower = 0
  )
  data <- data.table::as.data.table(data)
  data[, c("id", "time", "status") := list(
    seq_len(.N), max(cut), 1L
  )]

  form <- Formula::Formula(formula)
  f1 <- formula(form, rhs = 1)

  if (length(form)[[2]] > 1) {
    stop("Two-part formulas are not supported by sim_pexponential.")
  }

  atomic_cols <- colnames(data)[vapply(data, is.atomic, logical(1))]
  ped <- pammtools::split_data(
    formula = survival::Surv(time, status) ~ .,
    data = data[, ..atomic_cols],
    cut = cut,
    id = "id"
  )
  colnames(ped)[colnames(ped) == "tstart"] <- "t"

  eta <- lazyeval::f_eval(f1, data = ped)
  pammtools:::assert_numeric_eta(eta, f1, data)
  ped[["rate"]] <- exp(eta)
  ped <- data.table::as.data.table(ped)

  sim_data <- suppressMessages(ped[, .(
    time = pammtools:::rpexp(rate = rate, t = t),
    status = 1L
  ), by = id])
  sim_data <- data[, !c("time", "status")][sim_data, on = "id"]
  as.data.frame(sim_data)
}
