### setup
# Packages
source("setup.R")
# Functions
source("R/sim_pexponential.R")
source("R/sim_censoring.R")
source("R/simulation_pam.R")
source("R/pemtree_splines.R")

##### Main simulation

### Run info
run_id <- "compare_pemtree_cluster_argument_v1"
block_size <- 50

### Parameters
log_hazard <- function(t, delta) {
  -2 + delta * (t >= tau)
}
