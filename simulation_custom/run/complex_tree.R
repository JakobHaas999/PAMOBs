##########################################################################
## Complex tree DGP: can pemtree reconstruct a tree with 3 split levels? ##

### setup
source("setup.R")
source("R/sim_surv.R")
source("R/pemtree_splines.R")
library("parallel")

n_sim <- 10
n_obs <- 1000
n_pilot <- 5000

admin_time <- 10
censoring_rate <- 0.4

alpha <- 0.05
maxdepth <- 5

n_cores <- 6

## Baseline
shape <- 1.5
rate <- 0.05
log_baseline_hazard <- function(t) {
  log(shape) + log(rate) + (shape - 1) * log(t)
}

# Leaf assignments and leaf specified effects
true_leaf <- function(z1, z2, z3) {
  ifelse(
    z1 <= 0,
    ifelse(z2 <= 6, 1, ifelse(z3 <= 10, 2, 3)),
    ifelse(z3 <= 10, 4, 5)
  )
}
b1 <- c(0.5, 1, 1, -0.5, -0.5)
b2 <- c(0.5, -0.5, 0.5, 0, 1)

# sim_formula <-
