# Fit a ridge regression model to the data, and estimate the penalty
# strength (i.e., the normal prior on the regression coefficients)
# using cross-validation.
fit_ridge = function(X, y, lambda = c("lambda.min", "lambda.1se"),
                     standardize = FALSE)
  fit_lasso(X,y,lambda,standardize,0)

# Fit a Lasso model to the data, and estimate the penalty strength
# (lambda) using cross-validation.
fit_lasso = function(X, y, lambda = c("lambda.min", "lambda.1se"),
                     standardize = FALSE, alpha = 1) {
  timing = system.time(
    fit <- glmnet::cv.glmnet(x = X, y = y, standardize = standardize,
                             alpha = 1))
  b = drop(coef(fit,s = fit[[lambda]]))
  return(list(fit = fit,mu = b[1],beta = b[-1],timing = timing))
}

# Compute a fully-factorized variational approximation for Bayesian
# variable selection in linear regression.
fit_varbvs = function(X, y) {
  logodds = seq(-log10(ncol(X)),1,length.out = 40)
  timing = system.time(
    fit <- varbvs::varbvs(X, Z = NULL, y, logodds = logodds, verbose = FALSE))
  b = drop(coef(fit)[,"averaged"])
  return(list(fit = fit,mu = b[1],beta = b[-1],timing = timing))
}
