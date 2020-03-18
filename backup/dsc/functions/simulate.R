# Simulate training and test data from linear regression model y = Xb
# + e, where e is i.i.d. normal with zero mean and standard deviation
# "sigma". The outputs are the data (X, X.test), the outcomes (y,
# y.test), and the regression coefficients (beta) and residual
# s.d. (sigma) used to simulate the outcomes.
simulate_data = function(n = NULL, p = NULL, s = NULL, design = "indepgauss",
                         rho = 0, filepath = NULL, beta = NULL,
                         signal = "normal", sigma = NULL, Sigma.sqrt = NULL,
                         pve = 0.5, snr = NULL) {
  
  # Determine the simulation design.
  if (design == "indepgauss") {
    data              = list(X = matrix(rnorm(n*2*p), n*2, p))
  } else if (design == "equicorrgauss") {
    data              = list(X = matrix(rnorm(n*2*p), n*2, p))
    n                 = nrow(data$X)
    data$X            = sqrt(rho)*rnorm(n) + sqrt(1 - rho)*data$X
  } else if (design == "localcorrgauss") {
    data              = list(X = rmvnorm(2 * n, sigma = Sigma.sqrt))
  } else if (design == "realgenotype") {
    data              = readRDS(filepath)
    data              = list(X = scale(data$X))
    n                 = nrow(data$X) / 2
  } else if (design == "changepoint") {
    data              = list(X = matrix(0, n, n-1))
    for(j in 1:(n-1)) {
      for(i in (j+1):n) {
        data$X[i,j] = 1
      }
    }
  }
  
  # Simulate X and X.test.
  if (design == "changepoint") {
    X               = data$X
    X.test          = data$X
  } else {
    train.index     = sample(2*n, n)
    test.index      = setdiff(1:(2*n),train.index)
    X               = data$X[train.index,]
    X.test          = data$X[test.index,]
  }
  
  # Check whether X has "trivial" columns and, if so, remove them.
  var.X               = drop(apply(X, 2, var))
  ind                 = which(var.X / max(var.X) < 1/n)
  if (length(ind) > 0) {
    X                 = X[,-ind]
    X.test            = X.test[,-ind]
  }
  p                   = ncol(X)
  
  # Simulate the regression coefficients.
  if (is.null(beta)) {
    if (pve == 0) {
      beta          = double(p)
    } else {
      beta          = simulate_beta(p,s,signal = signal)
    }
  } else if (pve == 0)
    beta = double(p)  
  
  # Set the s.d. of the residual.
  if (is.null(sigma)) {
    if (pve == 0) {
      sigma          = 1
    } else {
      sigma          = sqrt(var(drop(X %*% beta)) * (1-pve) / pve)
    }
  }
  
  if (design == "changepoint") {
    sigma            = max(abs(beta)) / snr
  }
  
  # Sample the outcomes, y and y.test
  y                = X %*% beta + sigma * rnorm(n)
  err.test         = sigma * rnorm(n)
  if (design == "changepoint") {
    y.test         = X.test %*% beta
  } else {
    y.test         = X.test %*% beta + err.test
  }
  
  return (list(X = X, X.test = X.test, y = y, y.test = y.test,
               sigma = sigma, beta = beta))
}

# Simulate a (possibly sparse) vector of regression coefficients.
simulate_beta = function(p, s, s1 = 10, signal = "normal") {
  beta = double(p)
  if (signal == "t2") {
    beta[sample(p,s)] = rt(s, df = 2)
  } else if (signal == "t5") {
    beta[sample(p,s)] = rt(s, df = 5)
  } else if (signal == "lap") {
    beta[sample(p,s)] = rexp(s) * sign(rnorm(s))
  } else if (signal == "normal") {
    beta[sample(p,s)] = rnorm(s)
  } else if (signal == "unif") {
    beta[sample(p,s)] = runif(s)
  } else if (signal == "const") {
    beta[sample(p,s)] = 1
  } else if (signal == "subogdancandes") {
    ind               = sample(p,s)
    beta[ind]         = 0.1
    beta[ind[1:s1]]   = 1
  } else if (signal == "polygenic") {
    ind               = sample(p,s)
    beta[ind]         = 1
    beta[ind[1]]      = sqrt(s)
  }
  return (beta)
}

