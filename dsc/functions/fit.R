# Compute a fully-factorized variational approximation for Bayesian
# variable selection in linear regression. Input X should be an n x p
# numeric matrix, and input y should be a numeric vector of length n.
fit_varbvs = function(X, y) {
  logodds = seq(-log10(ncol(X)),1,length.out = 40)
  timing = system.time(
    fit <- varbvs::varbvs(X, Z = NULL, y, logodds = logodds, verbose = FALSE))
  b = drop(coef(fit)[,"averaged"])
  return(list(fit = fit,mu = b[1],beta = b[-1],timing = timing))
}
