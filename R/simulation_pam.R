### Function to simulate from a pam
#
# # Example call:
# n <- 1000
# covariate_spec <- list(
#   x1 = list(
#     distfun = rnorm
#   ),
#   x2 = list(
#     distfun = rbinom,
#     size = 1,
#     prob = 0.5
#   )
# )
# baseline <- function(t) {
#   ifelse(t <= 4, dgamma(t, 8, 2), log(0.5))
# }
# t_vals <- seq(0, 8, by = 0.01)
# plot(t_vals, exp(baseline(t_vals)),
#   type = "l",
#   xlab = "time",
#   ylab = expression(log(lambda[0](t))),
#   ylim = c(0, 2)
# )
# formula <- ~ baseline(t) + 0.5 * x1 + 0.8 * x2
# set.seed(2209)
# surv_data <- sim_pam(
#   n = n,
#   formula = formula,
#   covariate_spec = covariate_spec,
#   cut = seq(0, 8, by = 0.05)
# )
# ped <- pammtools::as_ped(
#   formula = Surv(time, status) ~ .,
#   data = surv_data,
#   cut = seq(0, 8, by = 0.1)
# )
# pam_fit <- mgcv::gam(
#   formula = ped_status ~ s(tend) + x1 + x2,
#   family = poisson(),
#   data = ped,
#   offset = ped$offset
# )
# coef(pam_fit)
# #> (Intercept)          x1          x2   s(tend).1   s(tend).2   s(tend).3
# #> 0.02097615  0.53853621  0.76870207 -0.56351083  5.26891387 -1.65825780
# #> s(tend).4   s(tend).5   s(tend).6   s(tend).7   s(tend).8   s(tend).9
# #> 0.73877876 -1.12479925  1.37489682  1.63042499 -7.08179795  1.11993399
# t_grid <- seq(0.1, 8, by = 0.1)
# newdata <- data.frame(
#   tend = t_grid, x1 = 0,  x2 = 0
# )
# hazard_hat <- predict(
#   pam_fit,
#   newdata = newdata,
#   type = "response"
# )
# hazard_true <- exp(baseline(t_grid))
# plot(t_grid, hazard_true, type = "l",
#      lwd = 2,
#      xlab = "t",
#      ylab = expression(lambda[0](t)),
#      ylim = range(hazard_hat, hazard_true))
# lines(t_grid, hazard_hat, lwd = 2, lty = 2)
# legend(
#   "topleft",
#   legend = c("true", "estimate"),
#   lwd = 2,
#   lty = c(1, 2)
# )

sim_pam <- function(n,
                    formula,
                    covariate_spec,
                    baseline,
                    censoring_rate,
                    cut) {
  # Arguments:
  # n: Number of observations for the simulation
  # formula:
  #   One-sided formula specifying the log-hazard.
  #   Functions of time and covariates are evaluated by
  #   pammtools::sim_pexp().
  # covariate_spec: List of Specifications of covariates
  # baseline: the log-hazard
  # censoring_rate: Proportion of censored observations
  # cut: A sequence of time-points starting with 0

  checkmate::assert_count(n)
  checkmate::assert_formula(formula)
  checkmate::assert_list(covariate_spec, null.ok = TRUE)
  checkmate::assert_function(baseline)
  checkmate::assert_numeric(cut, lower = 0, any.missing = FALSE)

  # Generate covariates
  if (!is.null(covariate_spec)) {
    covar <- lapply(covariate_spec, function(s) {
      checkmate::assert_list(s, names = "unique")
      checkmate::assert_true("distfun" %in% names(s))
      checkmate::assert_function(s$distfun)

      fun <- s$distfun
      args <- s[names(s) != "distfun"]

      do.call(fun, c(n = n, args))
    })

    data <- as.data.frame(covar)
    data$id <- seq_len(n)
  } else {
    data <- data.frame(id = seq_len(n))
  }

  sim_data <- pammtools::sim_pexp(
    formula = formula,
    data = data,
    cut = cut
  )

  censored_data <- sim_censoring(sim_data$time, censoring_rate = censoring_rate)

  cbind.data.frame(
    censored_data, data
  )
}
