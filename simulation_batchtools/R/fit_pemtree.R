fit_pemtree <- function(dat, pars) {
  do.call(pemtree, c(list(data = dat), pars$model_args))
}