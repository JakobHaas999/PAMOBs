### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
source("R/permtree.R")

## Helper functions ---------------------------------------------------
run_simulation <- function() {
                             ## Comparison of PErMtree and pemtree ---------------------------------
                             ## 3 scenarios are checked here
                             ## Scenario 1: Null scenario with no interactions
                             ## Scenario 2: Power Scenario with an effect of 0.5
                             ## Scenario 3: More complex scenario with two possible split variables
                             ##             -> One has an effect and one has no effect
  ## Global parameters ---------------------------------------------------
  n_pilot <- 10000
}
n <- 1000
shape <- 1.5
rate <- 0.05

log_hazrd_base <- function(t) {
  log(shape) + log(rate) + (shape - 1) * log(t)
}
formula_base <- ~ log_hazard_base(t)

admin_time <- 10
censoring_rate <- 0.3

n_sim <- 50

## Scenario 1 ---------------------------------------------------------
covariate_spec <- list(
  x = list(distfun = rnorm),
  z = list(distfun = rbinom, size = 1, prob = 0.5)
)
