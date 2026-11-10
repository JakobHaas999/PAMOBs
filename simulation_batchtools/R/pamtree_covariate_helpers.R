source("simulation_batchtools/R/pemtree_covariate_helpers.R")
source("simulation_batchtools/R/tree_comparison_helpers.R")

fit_and_summarise_pam <- function(instance, bs_df, gam_method, sp, cut_step, alpha, maxdepth, min_events) {
  data <- instance[, c("time", "status", "x", "z")]
  sp <- if (is.na(sp)) NULL else sp

  fits <- vector("list", length = 3L)
  fits[[1]] <- tryCatch(
    pemtree(
      formula = Surv(time, status) ~ bs(tend, df = bs_df) + x + offset(offset) | z,
      data = data,
      cut = seq(0, ceiling(max(data$time)), by = cut_step),
      min_events = min_events,
      alpha = alpha,
      maxdepth = maxdepth
    ),
    error = function(e) e
  )
  fits[[2]] <- tryCatch(
    pamtree(
      formula = Surv(time, status) ~ s(tend, bs = "ps", k = 10) + x + offset(offset) | z,
      data = data,
      cut = seq(0, ceiling(max(data$time)), by = cut_step),
      min_events = min_events,
      sp = sp,
      gam_method = gam_method,
      alpha = alpha,
      maxdepth = maxdepth
    ),
    error = function(e) e
  )
  fits[[3]] <- tryCatch(
    mob(
      formula = Surv(time, status) ~ x | z,
      data = data,
      fit = make_wbreg_fit(min_events = min_events),
      control = mob_control(alpha = alpha, maxdepth = maxdepth)
    ),
    error = function(e) e
  )

  summarise_fit <- function(fit, tree, model) {
    if (inherits(fit, "error")) {
      return(
        data.table(model = model, split = NA, p_root = NA_real_, split_point = NA_real_, error = conditionMessage(fit))
      )
    }

    root <- node_party(tree)
    split <- !is.terminal(root)
    p <- info_node(root)$p.value
    brk <- if (split) breaks_split(split_node(root)) else NULL
    data.table(
      model = model, split = split,
      p_root = if (length(p) == 1) p else NA_real_,
      split_point = if (is.null(brk)) NA_real_ else brk,
      error = NA_character_
    )
  }

  rbind(
    summarise_fit(fits[[1]], if (!inherits(fits[[1]], "error")) fits[[1]]$tree, "pemtree"),
    summarise_fit(fits[[2]], if (!inherits(fits[[2]], "error")) fits[[2]]$tree, "pamtree"),
    summarise_fit(fits[[3]], fits[[3]], "weibull_aft")
  )[, observed_censoring := mean(instance$status == 0)][]
}
