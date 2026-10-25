source("simulation_batchtools/designs/null_scenarios_base.R", local = TRUE)

make_level_shift_design <- function(
  n_rep,
  deltas = log(c(1.25, 1.5, 2, 3)),
  tau = 5
) {
  checkmate::assert_count(n_rep, positive = TRUE)
  checkmate::assert_numeric(
    deltas, min.len = 1L, any.missing = FALSE, finite = TRUE
  )
  checkmate::assert_number(tau, finite = TRUE, lower = 0)
  if (anyDuplicated(deltas)) stop("deltas must contain distinct values.")
  if (n_rep > 10000L) stop("The seed scheme supports at most 10000 repetitions.")

  base <- make_null_design(n_rep = 1L)
  if (tau <= 0 || any(tau >= base$admin_time)) {
    stop("tau must lie strictly between zero and admin_time.")
  }
  grid <- expand.grid(
    base_row = seq_len(nrow(base)),
    delta = deltas,
    KEEP.OUT.ATTRS = FALSE
  )

  settings <- base[grid$base_row, , drop = FALSE]
  settings$delta <- grid$delta
  settings$tau <- tau

  settings$setting_id <- seq_len(nrow(settings))

  design <- settings[
    rep(seq_len(nrow(settings)), each = n_rep), , drop = FALSE
  ]

  design$rep <- rep(seq_len(n_rep), times = nrow(settings))
  design$seed <- (
    2000000L + 10000L * design$setting_id + design$rep
  )
  rownames(design) <- NULL

  design
}
