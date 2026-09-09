source("setup.R")

# load data
tumor <- pammtools::tumor

tumor <- tumor[1:200, ]
tumor <- tumor %>%
  select(days, status, charlson_score, age, sex)

tumor_without_covar <- tumor %>% select(days, status)

tumor_ped <- as_ped(
  Surv(days, status) ~ 1,
  data = tumor_without_covar,
  cut = seq(0, max(tumor_without_covar$days), by = 0.5)
)

tree <- glmtree(
  formula = ped_status ~ 1 | tend,
  data = tumor_ped,
  offset = tumor_ped$offset,
  family = poisson(),
  alpha = 0.05,
  maxdepth = 1
)

plot(tree)
coef(tree)
