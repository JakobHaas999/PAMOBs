make_wbreg_fit <- function(min_events) {
  force(min_events)
  function(y, x, start = NULL, weights = NULL, offset = NULL, ..., estfun = FALSE, object = FALSE) {
    invalid <- list(
      coefficients = setNames(numeric(ncol(x)), colnames(x)),
      objfun = Inf, estfun = NULL, object = NULL
    )

    if (sum(y[, "status"]) < min_events) {
      return(invalid)
    }

    fit <- tryCatch(
      suppressWarnings(survreg(y ~ 0 + x, weights = weights, dist = "weibull")),
      error = function(e) NULL
    )
    if (is.null(fit) || any(!is.finite(coef(fit)))) {
      return(invalid)
    }

    list(
      coefficients = coef(fit),
      objfun = -as.numeric(logLik(fit)),
      estfun = if (estfun) estfun(fit) else NULL,
      object = if (object) fit else NULL
    )
  }
}

fit_and_summarise_tree_comp <- function(instance, bs_df, cut_step, alpha, maxdepth, min_events) {
  data <- instance[, c("time", "status", "x", "z")]

  fits <- vector("list", 2)
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
      return(data.table(
        model = model, split = NA, p_root = NA_real_, split_point = NA_real_,
        error = conditionMessage(fit)
      ))
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
    summarise_fit(fits[[2]], fits[[2]], "weibull_aft")
  )[, observed_censoring := mean(instance$status == 0)][]
}
