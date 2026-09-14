## source relevant files
source("setup.R")
source("R/simulation_covariates.R")
source("R/simulation_pem.R")
source("R/simulation_censoring.R")
source("R/simulation_pem_data.R")
source("R/fit_mob_pem.R")

## scenarios:
# 1. no interaction -> no split expected
# 2. interaction effect with a binary variables
# 3. interaction between two metric covariables

n <- 10000

## scenario 1
set.seed(1)
covariate_spec <- list(
  x1 = list(
    distfun = rnorm,
    mean = 0,
    sd = 1
  ),
  x2 = list(
    distfun = rbinom,
    size = 1,
    prob = 0.5
  )
)
formula <- ~ x1 + x2
lambdas <- c(0.1, 0.2, 0.3)
cuts <- c(2, 4)
beta <- c(0.8, 0)
censoring_rate <- 0.3

surv_data <- sim_pem_data(
  n = n,
  covariate_spec = covariate_spec,
  formula = formula,
  lambdas = lambdas,
  cuts = cuts,
  beta = beta,
  censoring_rate = censoring_rate
)

pem_formula <- Surv(time, status) ~ x1
z_formula <- ~x2
tree <- pemob(
  pem_formula = pem_formula,
  z_formula = z_formula,
  data = surv_data,
  cut = c(0, 2, 4, max(surv_data$time))
)
coef(tree)

## scenario 2
set.seed(2)
covariate_spec <- covariate_spec
formula <- ~ x1 + x2 + x1:x2
lambdas <- c(0.1, 0.2, 0.3)
cuts <- c(2, 4)
beta <- c(0.8, 0.5, 1.2)
censoring_rate <- 0.3

surv_data <- sim_pem_data(
  n = n,
  covariate_spec = covariate_spec,
  formula = formula,
  lambdas = lambdas,
  cuts = cuts,
  beta = beta,
  censoring_rate = censoring_rate
)

pem_formula <- Surv(time, status) ~ x1
z_formual <- ~x2

tree <- pemob(
  pem_formula = pem_formula,
  z_formula = z_formula,
  data = surv_data,
  cut = c(0, 2, 4, max(surv_data$time))
)
coef(tree)

## scenario 3
set.seed(3)
covariates <- sim_covariates(
  n = n,
  spec = list(
    x1 = list(
      distfun = rnorm,
      mean = 0,
      sd = 1
    ),
    x2 = list(
      distfun = runif,
      min = 0,
      max = 6
    )
  )
)
group <- as.integer(covariates$x2 > 5)
beta1 <- 0.3
beta2 <- 0.8

X <- cbind(
  x1 = covariates$x1,
  x1_group = covariates$x1 * group
)

beta <- c(beta1, beta2 - beta1)

times <- sim_pem(
  n = n,
  lambdas = c(0.1, 0.2, 0.3),
  cuts = c(2, 4),
  X = X,
  beta = beta
)

surv_data <- cbind(
  time = times,
  status = 1,
  covariates
)

tree <- pemob(
  pem_formula = Surv(time, status) ~ x1,
  z_formula = ~x2,
  data = surv_data,
  cut = c(0, 2, 4, max(surv_data$time))
)
plot(tree)
