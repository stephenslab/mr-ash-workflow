# sample possibly sparse beta
sample_beta = function(p, s, s1 = 10, signal = "normal") {
  
  beta              = double(p)
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

# simulate data
simulate_data = function(n = NULL, p = NULL, s = NULL, seed = 1,
                         design = "indepgauss", rho = 0,
                         filepath = NULL,
                         beta = NULL, signal = "normal",
                         sigma = NULL, Sigma.sqrt = NULL,
                         pve = 0.5, snr = NULL) {
  
  # set seed
  set.seed(2010 + seed)
  
  # design setting
  if (design == "indepgauss") {
    data              = list(X = matrix(rnorm(n*2*p), n*2, p))
  } else if (design == "equicorrgauss") {
    data              = list(X = matrix(rnorm(n*2*p), n*2, p))
    data$X           <- rnorm(dim(data$X)[1]) * sqrt(rho) + data$X * sqrt(1-rho)
  } else if (design == "identity") {
    data              = list(X = diag(n))
  } else if (design == "orthogonal") {
    cat("'2n' must be smaller than p\n")
    data              = list(X = qr.Q(qr(matrix(rnorm(n*2*p), n*2, p))))
  } else if (design == "localcorrgauss") {
    data              = list(X = rmvnorm(2 * n, sigma = Sigma.sqrt))
  } else if (design == "realgenotype") {
    data              = readRDS(filepath);
    data              = list(X = scale(data$X))
    n                 = dim(data$X)[1] / 2
  } else if (design == "changepoint") {
    data              = list(X = matrix(0, n, n-1))
    for(j in 1:(n-1)) {
      for(i in (j+1):n) {
        data$X[i,j] = 1
      }
    }
  }
  
  # sample X and X.test
  if ((design == "changepoint") || (design == "identity")) {
    X               = data$X
    X.test          = data$X
  } else {
    train.index       = sample(2*n, n)
    test.index        = (1:(2*n))[-train.index]
    X                 = data$X[train.index,]
    X.test            = data$X[test.index,]
  }
  
  # check whether X has trivial columns
  var.X               = c(apply(X, 2, var))
  ind                 = which(var.X / max(var.X) < 1 / n)
  if (length(ind) > 0) {
    X                 = X[,-ind]
    X.test            = X.test[,-ind]
  }
  p                   = dim(X)[2]
  
  # sample beta
  if (is.null(beta)) {
    if (pve == 0) {
      beta          = double(p)
    } else {
      beta          = sample_beta(p,s,signal = signal)
    }
  }
  
  # set sigma
  if (is.null(sigma)) {
    if (pve == 0) {
      sigma          = 1
    } else {
      sigma          = sqrt(var(c(X %*% beta)) * (1-pve) / pve)
    }
  }
  
  if (design == "changepoint") {
    sigma            = max(abs(beta)) / snr
  }
  
  # sample y and y.test
  y                <- X %*% beta + sigma * rnorm(n)
  err.test          = sigma * rnorm(n)
  if (design == "changepoint") {
    y.test           <- X.test %*% beta
  } else {
    y.test           <- X.test %*% beta + err.test
  }
  
  return (list(X = X, X.test = X.test, y = y, y.test = y.test,
               sigma = sigma, beta = beta))
}
