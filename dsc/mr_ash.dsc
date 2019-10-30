# A DSC for evaluating prediction accuracy of linear regression
# methods in different scenarios.
DSC:
  lib_path:  functions
  exec_path: modules/simulate,
             modules/fit,
             modules/predict,
             modules/score
  replicate: 20
  define:
    simulate: sparse_normal
  run: simulate

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
  $X:     out$train$X
  $y:     out$train$y
  $Xtest: out$X.test
  $ytest: out$y.test
  $beta:  out$beta
  $sigma: out$sigma
