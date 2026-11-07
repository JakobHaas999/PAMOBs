log_baseline <- function(t) {
  log(0.05) + log(1.5) + 0.5 * log(t)
}

generate_data <- function(n, scenario, beta, dz, dt, tau, admin_time, lambda_c) {
  x <- rbinom(n, size = 1, prob = 0.5)
  z <- round(rnorm(n, mean = 0, sd = 10))
  if (scenario == "covar_time_chang") {
    sim_formula <- ~ log_baseline(t) + beta * x + dz * x * (z > 0) + dt * x * (z > 0) * (t >= tau)
  } else if (scenario == "covar_struc_chang") {
    sim_formula <- ~ log_baseline(t) + beta * x + dz * x * (z > 0) + dt * (z > 0) * (t >= tau)
  } else {
    stop("Unknown scenario")
  }

  data <- data.frame(id = seq_len(n), x = x, z = z)

  sim_surv(
    sim_formula, data,
    breaks = tau, admin_time = admin_time, lambda_c = lambda_c
  )
}

# Plot oracle covariate interaction
# beta <- -1
# dz <- 1
# dt <- 1
# tau <- 3
#
# set.seed(2026)
# source("R/sim_surv.R")
# dat <- generate_data(20000, "covar_time_chang", beta, dz, dt, tau, admin_time = 10, lambda_c = 0.05)
#
# library("pammtools")
# library("survival")
# ped <- as_ped(Surv(time, status) ~ x + z, data = dat, cut = seq(0, 10, by = 0.1))
# library("data.table")
# setDT(ped)
# ped[, older := as.integer(z > 0)]
# ped[, post := as.integer(tend > tau + 1e-8)]
#
# oracle <- glm(ped_status ~ log(tend) + x + x:older + x:older:post + offset(offset),
#   data = ped, family = poisson()
# )
#
# est <- summary(oracle)$coefficients
# tab <- data.table(
#   term = rownames(est),
#   truth = c(log(0.05) + log(1.5), 0.5, beta, dz, dt),
#   estimate = est[, "Estimate"],
#   se = est[, "Std. Error"]
# )
# tab[, z_score := (estimate - truth) / se]
# tab
#
# grid <- CJ(tend = seq(0.1, 10, by = 0.1), x = 0:1, older = 0:1)
# grid[, post := as.integer(tend > tau + 1e-8)]
# grid[, offset := 0]
# grid[, oracle := predict(oracle, newdata = grid, type = "response")]
# grid[, truth := exp(log_baseline(tend) + x * (beta + older * (dz  + dt * post)))]
#
# library("ggplot2")
# ggplot(grid, aes(x = tend, colour = factor(x))) +
#   geom_line(aes(y = oracle, linetype = "Oracle fit"), linewidth = 0.7) +
#   geom_line(aes(y = truth, linetype = "Truth"), linewidth = 0.7) +
#   geom_vline(xintercept = tau, colour = "grey75", linewidth = 0.4) +
#   facet_wrap(~ older, labeller = labeller(older = c("0" = "z <= 0", "1" = "z > 0"))) +
#   scale_colour_manual(values = c("0" = "#0072B2", "1" = "#D55E00"), labels = c("Control", "Treated"), name = NULL) +
#   scale_linetype_manual(values = c("Oracle fit" = "solid", "Truth" = "22"), name = NULL) +
#   labs(x = "Time", y = expression(lambda(t))) +
#   theme_bw()

calibrate_one <- function(n_calibration, scenario, beta, dz, dt, tau, admin_time, censoring_rate) {
  d <- generate_data(
    n = n_calibration,
    scenario = scenario,
    beta = beta,
    dz = dz,
    dt = dt,
    tau = tau,
    admin_time = admin_time,
    lambda_c = 0
  )
  event_times <- data.frame(eventtime = d$time)

  data.table(
    censoring_rate = censoring_rate,
    lambda_c = vapply(censoring_rate, function(rate) {
      calibrate_censoring(event_times, rate, admin_time)
    }, numeric(1))
  )
}

get_splits <- function(node, vars, path = "root") {
  if (is.terminal(node)) {
    return(NULL)
  }
  s <- split_node(node)
  p <- info_node(node)$p.value
  rbind(
    data.table(
      path = path, var = vars[varid_split(s)], at = breaks_split(s),
      p = if (length(p) == 1) unname(p) else NA_real_
    ),
    get_splits(kids_node(node)[[1]], vars, paste0(path, "L")),
    get_splits(kids_node(node)[[2]], vars, paste0(path, "R"))
  )
}

fit_and_summarise_complex_tree <- function(instance, cut_step, alpha, maxdepth, min_events) {
  data <- instance[, c("time", "status", "x", "z")]

  fit_args <- list(
    formula = Surv(time, status) ~ bs(tend, df = 6) + x + offset(offset) | tend + z,
    data = data,
    min_events = min_events,
    alpha = alpha,
    maxdepth = maxdepth
  )

  if (!is.na(cut_step)) {
    fit_args$cut <- seq(0, ceiling(max(data$time)), by = cut_step)
  }
  fit <- tryCatch(do.call(pemtree, fit_args), error = function(e) e)

  sp <- data.table(path = character(), var = character(), at = numeric(), p = numeric())
  error <- NA_character_
  p_root <- NA_real_

  if (inherits(fit, "error")) {
    error <- conditionMessage(fit)
  } else {
    tree <- fit$tree
    root <- node_party(tree)
    p <- info_node(root)$p.value
    if (length(p) == 1) p_root <- unname(p)
    sp_fit <- get_splits(root, colnames(tree$data))
    if (!is.null(sp_fit)) sp <- sp_fit
  }

  node_split <- function(where) {
    s <- sp[path == where]
    if (nrow(s) == 0) {
      return(list(var = "none", at = NA_real_, p = NA_real_))
    }
    list(var = s$var, at = s$at, p = s$p)
  }

  root_s <- node_split("root")
  right_s <- node_split("rootR")
  left_s <- node_split("rootL")

  data.table(
    n_split = if (is.na(error)) nrow(sp) else NA_integer_,
    root_var = root_s$var, root_at = root_s$at, p_root = p_root,
    right_var = right_s$var, right_at = right_s$at, p_right = right_s$p,
    left_var = left_s$var, left_at = left_s$at, p_left = left_s$p,
    splits = list(sp),
    error = error,
    observed_censoring = mean(instance$status == 0)
  )
}
