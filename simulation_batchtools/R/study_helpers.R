# Cross data-generating settings with model df and repetitions. Seeds identify
# data sets, so the same data are regenerated for every df in a comparison.
cross_df_repetitions <- function(settings, dfs, n_rep, seed_offset) {
  settings$dgp_id <- seq_len(nrow(settings))
  grid <- expand.grid(data_row = seq_len(nrow(settings)), df = as.integer(dfs))
  settings <- settings[grid$data_row, , drop = FALSE]
  settings$df <- grid$df
  settings$setting_id <- seq_len(nrow(settings))
  design <- settings[rep(seq_len(nrow(settings)), each = n_rep), , drop = FALSE]
  design$rep <- rep(seq_len(n_rep), times = nrow(settings))
  design$seed <- as.integer(seed_offset + 10000 * design$dgp_id + design$rep)
  rownames(design) <- NULL
  design
}

prepare_spline_pars <- function(pars) {
  model_terms <- c(
    sprintf("splines::bs(tend, df = %d)", pars$df),
    if (pars$covariate_setting == "included") c("x1", "x2"),
    "offset(offset)"
  )
  model_formula <- stats::as.formula(
    paste("survival::Surv(time, status) ~", paste(model_terms, collapse = " + "), "| tend"),
    env = asNamespace("stats")
  )
  pars$model_args <- list(
    formula = model_formula,
    cut = sort(unique(c(seq(0, pars$admin_time, by = pars$cut_step), pars$admin_time))),
    min_events = pars$min_events, alpha = pars$alpha, maxdepth = pars$maxdepth
  )
  pars
}
