### Simulate censoring times with a given rate
#
## Example function call
# set.seed(123)
# times <- sim_pem(
#   n = 1000,
#   lambdas = c(0.1, 0.3, 0.05),
#   cuts = c(2, 5)
# )
# censoring_rate <- 0.6
# surv_data <- sim_censoring(times, censoring_rate)
# head(surv_data)
# #> time status
# #> 1 1.0250338      0
# #> 2 0.5394712      0
# #> 3 9.5810974      1
# #> 4 0.3157736      1
# #> 5 0.5621098      1
# #> 6 1.4354894      0
# mean(surv_data$status == 0)
# #> [1] 0.604

sim_censoring <- function(times, censoring_rate) {
  # Arguments:
  # times: numeric vector of event times
  # censoring_rate: Proportion of censored observations

  # asserts
  checkmate::assert_numeric(times, any.missing = FALSE, lower = 0)
  checkmate::assert_true(all(t > 0))
  checkmate::assert_number(censoring_rate)
  checkmate::assert_true(censoring_rate > 0 && censoring_rate < 1)

  lambda_c <- uniroot(
    f = function(lambda_c) {
      p_censored <- mean(1 - exp(-lambda_c * times))
      p_censored - censoring_rate
    },
    interval = c(1e-8, 100)
  )$root

  censoring_times <- rexp(length(times), rate = lambda_c)

  data.frame(
    time = pmin(times, censoring_times),
    status = as.integer(times <= censoring_times)
  )
}
