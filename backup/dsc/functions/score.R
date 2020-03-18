# Return the square root of the sum of squared errors between the observed
# values ("obs") and the predicted, or "fitted," values.
rsse <- function (obs, fitted)
  norm((obs - fitted),'2')
