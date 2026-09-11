### Simulate covariates for a PEM or PAM
#
## Example call
# set.seed(123)
# n <- 1000
# spec <- list(
#   age = list(
#     distfun = rnorm,
#     mean = 44,
#     sd = 8
#   ),
#   treatment = list(
#     distfun = rbinom,
#     size = 1,
#     prob = 0.5
#   )
# )
# x <- sim_covariates(n, spec)
# head(x)
#> age treatment
#> 1 39.51619         0
#> 2 42.15858         0
#> 3 56.46967         0
#> 4 44.56407         1
#> 5 45.03430         0
#> 6 57.72052         1
sim_covariates <- function(n, spec) {
  # Arguments:
  # n: Number of simulations/observations
  # spec: List of Specifications of covariates



  checkmate::assert_count(n, positive = TRUE)
  checkmate::assert_list(spec, min.len = 1)

  covariates <- lapply(spec, function(s) {
    checkmate::assert_list(s, names = "unique")
    checkmate::assert_true("distfun" %in% names(s))

    distfun <- s$distfun
    args <- s[names(s) != "distfun"]

    do.call(distfun, c(list(n = n), args))
  })

  as.data.frame(covariates)
}
n <- 1000
