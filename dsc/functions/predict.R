# This file contains functions to predict outcomes given various
# fitted models.

# Predict outcomes y from a simple linear regression, y = mu + X*beta,
# where X is an n x p matrix of observations, beta is a vector of
# regression coefficients of length p, and mu is the intercept.
predict_linear <- function (X, mu, beta)
  drop(mu + X %*% beta)
