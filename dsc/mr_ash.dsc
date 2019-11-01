# A DSC for evaluating prediction accuracy of linear regression
# methods in different scenarios.
DSC:
  R_libs:    glmnet, varbvs >= 2.6-3
  lib_path:  functions
  exec_path: modules/simulate,
             modules/fit,
             modules/predict,
             modules/score
  replicate: 4
  define:
    simulate: sparse_normal, sparse_t
    fit:      ridge, lasso, varbvs
    predict:  predict_linear
    score:    rsse
  run: simulate * fit * predict * score

# simulate modules
# ================
# A "simulate" module generates a training and test data set used to
# evaluate each of the linear regression models. Each training and
# test data set should include an n x p matrix X and a vector y of
# length n, where n is the number of samples, and p is the number of
# candidate predictors.

# Simulate the linear regression coefficients such that the s nonzero
# coefficients are drawn from the normally distributed with zero
# mean. The residual variance is controlled to achieve the target PVE.
sparse_normal: sparse_normal.R
  n:      500
  p:      2000
  s:      1, 5, 20, 100, 500, 2000
  pve:    0.5
  $X:     out$X
  $y:     out$y
  $Xtest: out$X.test
  $ytest: out$y.test
  $beta:  out$beta
  $sigma: out$sigma

# Simulate the linear regression coefficients such that the s nonzero
# coefficients are drawn from the t distrubution with 2 degrees of
# freedom. The residual variance is controlled to achieve the target
# PVE.
sparse_t: sparse_t.R
  n:      500
  p:      2000
  s:      1, 5, 20, 100, 500, 2000
  pve:    0.5
  $X:     out$X
  $y:     out$y
  $Xtest: out$X.test
  $ytest: out$y.test
  $beta:  out$beta
  $sigma: out$sigma

# fit modules
# ===========
# A "fit" module fits a linear regression model to the provided
# training data, X and y.

# Fit a ridge regression model using glmnet. The penalty strength
# (i.e., the normal prior on the coefficients) is estimated using
# cross-validation.
ridge: ridge.R
  standardize:       FALSE
  lambda_est_method: "lambda.min", "lambda.1se"
  X:                 $X
  y:                 $y
  $intercept:        out$mu
  $beta_est:         out$beta
  $timing:           out$timing
  $model:            out$fit

# Fit a Lasso model using glmnet. The penalty strength ("lambda") is
# estimated via cross-validation.
lasso: lasso.R
  standardize:       FALSE
  lambda_est_method: "lambda.min", "lambda.1se"
  X:                 $X
  y:                 $y
  $intercept:        out$mu
  $beta_est:         out$beta
  $timing:           out$timing
  $model:            out$fit
  
# Compute a fully-factorized variational approximation for Bayesian
# variable selection in linear regression ("varbvs").
varbvs: varbvs.R
  X:          $X
  y:          $y
  $intercept: out$mu
  $beta_est:  out$beta
  $timing:    out$timing
  $model:     out$fit

# predict modules
# ===============
# A "predict" module takes as input a fitted model (or the parameters
# of this fitted model) and an n x p matrix of observations, X, and
# returns a vector of length n containing the outcomes predicted by
# the fitted model.

# Predict outcomes from a fitted linear regression model.
predict_linear: predict_linear.R
  X:         $Xtest
  intercept: $intercept
  beta:      $beta_est
  $yest:     y

# score modules
# =============
# A "score" module takes as input a vector of predicted outcomes and a
# vector of true outcomes, and outputs a summary statistic that can be
# used to evaluate accuracy of the predictions.

# Compute the square root of the sum of squared errors summarizing the
# differences between the predicted outcomes and the true outcomes.
rsse: rsse.R
  y:    $ytest
  yest: $yest
  $err: err
