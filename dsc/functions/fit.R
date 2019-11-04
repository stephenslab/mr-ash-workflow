# Fit a ridge regression model to the data, and estimate the penalty
# strength (i.e., the normal prior on the regression coefficients)
# using cross-validation.
fit_ridge = function(X, y, lambda = c("lambda.min", "lambda.1se"),
                     standardize = FALSE)
  fit_glmnet(X, y, match.arg(lambda), standardize, alpha = 0)

# Fit a Lasso model to the data (or, more generally, an Elastic Net
# model with a fixed setting of alpha), and estimate the penalty
# strength (lambda) using cross-validation.
fit_lasso = function(X, y, lambda = c("lambda.min", "lambda.1se"),
                     standardize = FALSE)
  fit_glmnet(X, y, match.arg(lambda), standardize, alpha = 1)

# Fit an Elastic Net model to the data, and estimate the Elastic Net
# parameters (penalty strength, "lambda", and mixing parameter,
# "alpha") using cross-validation.
fit_elastic_net = function(X, y, lambda = c("lambda.min", "lambda.1se"),
                           standardize = FALSE, alpha = seq(0,1,0.1)) {
  lambda = match.arg(lambda)

  # Fit an Elastic Net model for each candidate setting of the mixing
  # parameter, alpha.
  n      = length(alpha)
  out    = vector("list",n)
  cv.err = double(n)
  timing = system.time(
    for (i in 1:n) {
      out[[i]]   = fit_glmnet(X, y, lambda, standardize, alpha[i])
      fit        = out[[i]]$fit
      lambda.ind = which(fit$lambda == fit$lambda.1se)
      cv.err[i]  = fit$cvm[lambda.ind]
  })

  # Return the results for the elastic net model producing the lowest 
  i = which.min(cv.err)
  return(c(out[[i]][c("fit","mu","beta")],
           list(timing = timing["elapsed"],
                alpha  = alpha[i])))
}

# Fit an Elastic Net model, in which the penalty strength (lambda) is
# estimated using cross-validation.
fit_glmnet = function(X, y, lambda, standardize, alpha) {
  timing = system.time(
    fit <- glmnet::cv.glmnet(x = X, y = y, standardize = standardize,
                             alpha = alpha))
  b = drop(coef(fit, s = fit[[lambda]]))
  return(list(fit = fit, mu = b[1], beta = b[-1], timing = timing["elapsed"]))
}

# Compute a fully-factorized variational approximation for Bayesian
# variable selection in linear regression.
fit_varbvs = function(X, y, logodds = seq(-log10(ncol(X)),1,length.out = 40)) {
  timing = system.time(
    fit <- varbvs::varbvs(X, Z = NULL, y, logodds = logodds, verbose = FALSE))
  b = drop(coef(fit)[,"averaged"])
  return(list(fit = fit, mu = b[1], beta = b[-1], timing = timing["elapsed"]))
}

# Fit a fully-factorized variational approxition for Bayesian variable
# selection in linear regression, with scale-mixture-of-normals
# ("adaptive shrinkage") priors on the regression coefficients.
fit_mr_ash = function(X, y, standardize = FALSE, sa2 = (2^((0:19)/5) - 1)^2) {
  timing = system.time(
    fit <- mr.ash.alpha::mr.ash(X = X, y = y, sa2 = sa2, max.iter = 2000,
                                tol = list(epstol = 1e-12, convtol = 1e-8),
                                standardize = standardize))
  b = coef(fit)
  return(list(fit = fit, mu = b[1], beta = b[-1], timing = timing["elasped"]))
}
