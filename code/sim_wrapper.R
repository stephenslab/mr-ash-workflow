#'
#'
#'
#' sample possibly sparse beta
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
  } else if (signal == "const1") {
    ind               = sample(p,s)
    beta[ind]         = 0.1
    beta[ind[1:s1]]   = 1
  } else if (signal == "subogdancandes") {
    ind               = sample(p,s)
    beta[ind]         = 0.1
    beta[ind[1:s1]]   = 1
  } else if (signal == "polygenic") {
    ind               = sample(p,s)
    beta[ind]         = sqrt(p-s)
    beta[-ind]        = sqrt(s)
  }
  
  return (beta)
}

#'
#'
#'
#' simulate data
simulate_data = function(n = NULL, p = NULL, s = NULL, seed = 1,
                         design = "indepgauss",
                         filepath = NULL,
                         beta = NULL, signal = "normal",
                         sigma = NULL, pve = 0.5, rho = 0) {
  
  # set seed
  set.seed(2010 + seed)
  
  # design setting
  if (design == "indepgauss") {
    data              = list(X = matrix(rnorm(n*2*p), n*2, p))
  } else if (design == "equicorrgauss") {
    data              = list(X = matrix(rnorm(n*2*p), n*2, p))
    data$X           <- rnorm(dim(data$X)[1]) * sqrt(rho) + data$X * sqrt(1-rho)
  } else if (design == "realgenotype") {
    data              = readRDS(filepath);
    data              = list(X = scale(data$X))
    n                 = dim(data$X)[1] / 2
    p                 = dim(data$X)[2]
  } else if (design == "changepoint") {
    data$X            = list(X = zeros(0, n, n-1))
    for(j in 1:(n-1)) {
      for(i in (j+1):n) {
        data$X[i,j] = 1
      }
    }
  }
  
  # sample X and X.test
  if (design == "changepoint") {
    X               = data$X
    X.test          = data$X
  } else {
    train.index       = sample(2*n, n)
    test.index        = (1:(2*n))[-train.index]
    X                 = data$X[train.index,]
    X.test            = data$X[test.index,]
  }
  
  # sample beta
  if (is.null(beta)) {
    beta            = sample_beta(p,s,signal = signal)
  }
  
  # set sigma
  if (is.null(sigma)) {
    sigma             = sqrt(var(c(X %*% beta)) * (1-pve) / pve)
  }
  
  # sample y and y.test
  y                <- X %*% beta + sigma * rnorm(n)
  err.test          = sigma * rnorm(n)
  y.test           <- X.test %*% beta + err.test
  
  return (list(X = X, X.test = X.test, y = y, y.test = y.test,
               sigma = sigma, beta = beta))
}

#'
#'
#'
#' Plotting routines

## function for boxplot
my.box <- function (dat, x, y,
                    values = c(1,2,0,3,4,5)) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           scale_shape_manual(values = values) +
           geom_boxplot(alpha = 0.1, aes(color = fit), outlier.alpha = 0) +
           scale_alpha_manual(values = 0.1) +
           scale_fill_discrete(guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 14))
}

## function for line
my.line <- function (dat, x, y, cols = NULL,
                     shapes = NULL) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           scale_shape_manual(values = shapes) +
           scale_color_manual(values = cols) +
           geom_point(aes(color = fit, shape = fit)) +
           geom_line(aes(color = fit)) +
           scale_fill_discrete(guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 14))
}

## function for violin plot
my.jitter <- function (dat, x, y,
                       values = c(1,2,0,3,4,5)) {
  return(ggplot(dat,aes_string(x = x, y = y, fill = "fit")) +
           geom_jitter(aes(color = fit, shape = fit)) +
           scale_shape_manual(values = values) +
           geom_violin(alpha = 0.1, aes(color = fit, fill = fit), scale = "width") +
           scale_alpha_manual(0.1) +
           scale_fill_discrete(guide = "none") +
           labs(x           = "") +
           theme_cowplot(font_size = 14))
}

## function for color
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}