# A DSC for evaluating prediction accuracy of linear regression
# methods in different scenarios.
DSC:
  R_libs:    MASS, glmnet, susieR, varbvs >= 2.6-3
  lib_path:  functions
  exec_path: modules/simulate,
             modules/fit,
             modules/predict,
             modules/score
  replicate: 200
  define:
    simulate: null_effects, one_effect, zh
    fit:      ridge, lasso, elastic_net, susie, varbvs, varbvsmix
    predict:  predict_linear
    score:    mse, mae
  run: simulate * fit * predict * score

# simulate modules
# ================
# A "simulate" module generates a training and test data set used to
# evaluate the linear regression models. Each training and test data
# set should include an n x p matrix X and a vector y of length n,
# where n is the number of samples, and p is the number of candidate
# predictors.

# Simulate data in the same way as Example 1 of Zou & Hastie (2005),
# except that all the regression coefficients are zero.
null_effects: null_effects.R
  $X:     X
  $y:     y
  $Xtest: Xtest
  $ytest: ytest
  $beta:  b
  $se:    se

# Simulate data in the same way as Example 1 of Zou & Hastie (2005),
# except that all the regression coefficients except one are zero.
# The inputs, outputs and module parameters are the same as the "null"
# simulate module.
one_effect: one_effect.R
  $X:     X
  $y:     y
  $Xtest: Xtest
  $ytest: ytest
  $beta:  b
  $se:    se
  
# Generate training and test data sets using one of the four scenarios
# described in Zou & Hastie (2005). This module inherits the inputs,
# outputs and module parameters from the "null" simulate module.
zh: zh.R
  scenario: 1, 2, 3, 4
  $X:       out$train$X
  $y:       out$train$y
  $Xtest:   out$test$X
  $ytest:   out$test$y
  $beta:    out$b
  $se:      out$se

# fit modules
# ===========
# A "fit" module fits a linear regression model to the provided
# training data, X and y.

# Fit a ridge regression model using glmnet. The penalty strength
# (i.e., the normal prior on the coefficients) is estimated using
# cross-validation.
ridge: ridge.R
  X:          $X
  y:          $y
  $intercept: out$mu
  $beta_est:  out$beta
  $model:     out
  
# Fit a Lasso model using glmnet. The penalty strength ("lambda") is
# estimated via cross-validation.
lasso: lasso.R
  X:          $X
  y:          $y
  $intercept: out$mu
  $beta_est:  out$beta
  $model:     out

# Fit an Elastic Net model using glmnet. The model parameters, lambda
# and alpha, are estimated using cross-validation.
elastic_net: elastic_net.R
  X:          $X
  y:          $y
  $intercept: out$mu
  $beta_est:  out$beta
  $model:     out
    
# Fit a "sum of single effects" (SuSiE) regression model.
susie: susie.R
  X:          $X
  y:          $y
  $intercept: out$mu
  $beta_est:  out$beta
  $model:     out$fit

# Compute a fully-factorized variational approximation for Bayesian
# variable selection in linear regression (varbvs).
varbvs: varbvs.R
  X:          $X
  y:          $y
  $intercept: out$mu
  $beta_est:  out$beta
  $model:     out$fit

# This is a variant on the varbvs method in which the "spike-and-slab"
# prior on the regression coefficients is replaced with a
# mixture-of-normals prior.
varbvsmix: varbvsmix.R
  X:          $X
  y:          $y
  $intercept: out$mu
  $beta_est:  out$beta
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

# Compute the mean squared error summarizing the differences between
# the predicted outcomes and the true outcomes.
mse: mse.R
  y:    $ytest
  yest: $yest
  $err: err

# Compute the mean absolute error summarizing the differences between
# the predicted outcomes and the true outcomes.
mae: mae.R
  y:    $ytest
  yest: $yest
  $err: err
