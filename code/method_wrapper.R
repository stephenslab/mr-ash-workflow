#'
#'
#'
#' MR.ASH
fit.mr.ash = function(X, y, X.test, y.test, seed = 1, sa2 = (2^((0:19)/5) - 1)^2) {
  
  # set seed
  set.seed(seed)
  
  # run mr.ash
  t.mr.ash           = system.time(
    fit.mr.ash        <- mr.ash(X = X, y = y, sa2 = sa2,
                                max.iter = 2000,
                                standardize = standardize,
                                tol = list(epstol = 1e-12, convtol = 1e-8)))
  
  return (list(fit = fit.mr.ash, t = t.mr.ash[3],
               rsse = norm(y.test - predict(fit.mr.ash, X.test), '2')))
}

#'
#'
#'
#' MR.ASH
fit.mr.ash2 = function(X, y, X.test, y.test, seed = 1,
                       sa2 = (2^((0:19)/5) - 1)^2,
                       update.pi = TRUE, pi = NULL,
                       beta.init = NULL, sigma2 = NULL,
                       update.order = NULL) {
  
  
  # set seed
  set.seed(seed)
  
  # run mr.ash order
  t.mr.ash           = system.time(
    fit.mr.ash        <- mr.ash(X = X, y = y, sa2 = sa2, update.order = update.order,
                                max.iter = 2000, min.iter = 200,
                                beta.init = beta.init, update.pi = update.pi, pi = pi,
                                standardize = standardize, sigma2 = sigma2,
                                tol = list(epstol = 1e-12, convtol = 1e-8)))
  
  return (list(fit = fit.mr.ash, t = t.mr.ash[3],
               rsse = norm(y.test - predict(fit.mr.ash, X.test), '2')))
}

#'
#'
#'
#' Lasso
fit.lasso = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run lasso
  t.lasso           = system.time(
    fit.lasso        <- cv.glmnet(x = X, y = y, standardize = standardize))
  
  t.lasso2          = system.time(glmnet(x = X, y = y, standardize = standardize))
  
  return (list(fit = fit.lasso, t = t.lasso[3], t2 = t.lasso2[3],
               rsse = norm(y.test - predict(fit.lasso, newx = X.test, s = fit.lasso$lambda.1se), '2'),
               rsse2 = norm(y.test - predict(fit.lasso, newx = X.test, s = fit.lasso$lambda.min), '2')))
}


#'
#'
#'
#' Ridge
fit.ridge = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run lasso
  t.ridge           = system.time(
    fit.ridge        <- cv.glmnet(x = X, y = y, alpha = 0, standardize = standardize))
  
  return (list(fit = fit.ridge, t = t.ridge[3],
               rsse = norm(y.test - predict(fit.ridge, newx = X.test, s = fit.ridge$lambda.1se), '2'),
               rsse2 = norm(y.test - predict(fit.ridge, newx = X.test, s = fit.ridge$lambda.min), '2')))
}

fit.ridge.opt = function(X, y, X.test, y.test, sigma, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run ridge
  fit          = glmnet(x = X, y = (y - mean(y)) / sqrt(mean((y - mean(y))^2)),
                        alpha = 0, lambda = sigma^2 / length(y), standardize = standardize)
  fit$a0       = fit$a0 * sqrt(mean((y - mean(y))^2)) + mean(y)
  fit$beta     = fit$beta * sqrt(mean((y - mean(y))^2))
  
  return (list(fit = fit, t = 0,
               rsse = norm(y.test - predict(fit, newx = X.test, s = fit$lambda.1se), '2')))
}

#'
#'
#'
#' Elastic Net
fit.enet = function(X, y, X.test, y.test, seed = 1, all = FALSE) {
  
  # set seed
  set.seed(seed)
  
  # store CV error
  cv.err            = double(11)
  rsse              = double(11)
  rsse2             = double(11)
  time              = double(11)
  fit.enet          = list()
  
  # run ENET
  for (j in 0:10) {
    alpha           = 0.1 * j;
    t.enet          = system.time(
      fit.enet[[j+1]] <- cv.glmnet(x = X, y = y, alpha = alpha, standardize = standardize))
    lambda.ind      = which(fit.enet[[j+1]]$lambda == fit.enet[[j+1]]$lambda.1se)
    cv.err[j+1]     = fit.enet[[j+1]]$cvm[lambda.ind]
    rsse[j+1]       = norm(y.test - predict(fit.enet[[j+1]], newx = X.test,
                                            s = fit.enet[[j+1]]$lambda.1se), '2')
    rsse2[j+1]      = norm(y.test - predict(fit.enet[[j+1]], newx = X.test,
                                            s = fit.enet[[j+1]]$lambda.min), '2')
    time[j+1]       = as.numeric(t.enet[3])
  }
  
  alpha.ind         = which.min(cv.err)
  
  if (all == TRUE) {
    return (list(fit = fit.enet, t = time, rsse = rsse))
  }
  
  return (list(fit = fit.enet[[alpha.ind]], t = sum(time),
               rsse = rsse[alpha.ind],
               rsse2 = rsse2[alpha.ind],
               alpha = (alpha.ind - 1) * 0.1))
}

#'
#'
#'
#' MCP
fit.scad = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run SCAD
  gamma            = seq(2.1, 5.3, length.out = 11)
  for (j in 1:11) {
    t.scad2          = system.time(
      fit.scad2       <- cv.ncvreg(X, y, penalty = "SCAD", gamma = gamma[j], nfolds = 10))
    if (j == 1) {
      t.scad        <- t.scad2
      fit.scad      <- fit.scad2
    } else {
      t.scad        <- t.scad + t.scad2
      if (min(fit.scad$cve) > min(fit.scad2$cve)) {
        fit.scad      <- fit.scad2
      }
    }
  }
  
  return (list(fit = fit.scad, t = t.scad[3], rsse = norm(y.test - predict(fit.scad, X.test), '2')))
}

#'
#'
#'
#' MCP
fit.mcp = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run MCP
  gamma            = seq(1.1, 4.9, length.out = 11)
  for (j in 1:11) {
    t.mcp2           = system.time(
      fit.mcp2        <- cv.ncvreg(X, y, penalty = "MCP", gamma = gamma[j], nfolds = 10))
    if (j == 1) {
      t.mcp         <- t.mcp2
      fit.mcp       <- fit.mcp2
    } else {
      t.mcp         <- t.mcp + t.mcp2
      if (min(fit.mcp$cve) > min(fit.mcp2$cve)) {
        fit.mcp       <- fit.mcp2
      }
    }
  }
  
  return (list(fit = fit.mcp, t = t.mcp[3], rsse = norm(y.test - predict(fit.mcp, X.test), '2')))
}

#'
#'
#'
#' SCAD with a fixed gamma
fit.scad2 = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run SCAD
  t.scad            = system.time(
    fit.scad       <- cv.ncvreg(X, y, penalty = "SCAD", nfolds = 10))
  
  t.scad2           = system.time(ncvreg(X, y, penalty = "SCAD"))
  
  return (list(fit = fit.scad, t = t.scad[3], rsse = norm(y.test - predict(fit.scad, X.test), '2'),
               t2 = t.scad2[3]))
}

#'
#'
#'
#' MCP with a fixed gamma
fit.mcp2 = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run MCP
  t.mcp           = system.time(
    fit.mcp      <- cv.ncvreg(X, y, penalty = "MCP", nfolds = 10))
  t.mcp2          = system.time(ncvreg(X, y, penalty = "MCP"))
  
  return (list(fit = fit.mcp, t = t.mcp[3], rsse = norm(y.test - predict(fit.mcp, X.test), '2'),
               t2 = t.mcp2[3]))
}


#'
#'
#'
#' L0Learn
fit.l0learn = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run l0learn
  t.l0learn         = system.time(
    fit.l0learn      <- L0Learn.cvfit(X, y, nFolds = 10))
  lambda.min        = fit.l0learn$fit$lambda[[1]][which.min(fit.l0learn$cvMeans[[1]])]
  
  return (list(fit = fit.l0learn, t = t.l0learn[3],
               rsse = norm(y.test - predict(fit.l0learn, X.test, lambda = lambda.min)@x, '2')))
}



#'
#'
#'
#' Bayesian Lasso
fit.blasso = function(X, y, X.test, y.test, seed = 1, nIter = NULL, burnIn = NULL) {
  
  # set seed
  set.seed(seed)
  
  # default setting
  if (is.null(nIter)) nIter = 1500
  if (is.null(burnIn)) burnIn = 500
  
  # run blasso
  t.blasso          = system.time(
    fit.blasso       <- BGLR(y, ETA = list(list(X = X, model="BL", standardize = standardize)),
                             verbose = FALSE, nIter = nIter, burnIn = burnIn))
  fit.blasso$beta   = c(fit.blasso$ETA[[1]]$b)
  
  return (list(fit = fit.blasso, t = t.blasso[3],
               rsse = norm(y.test - X.test %*% fit.blasso$beta - fit.blasso$mu, '2')))
}


#'
#'
#'
#' BayesB
fit.bayesb = function(X, y, X.test, y.test, seed = 1, nIter = NULL, burnIn = NULL) {
  
  # set seed
  set.seed(seed)
  
  # default setting
  if (is.null(nIter)) nIter = 1500
  if (is.null(burnIn)) burnIn = 500
  
  # run bayesb
  t.bayesb          = system.time(
    fit.bayesb       <- BGLR(y, ETA = list(list(X = X, model="BayesB", standardize = standardize)),
                             verbose = FALSE, nIter = nIter, burnIn = burnIn))
  fit.bayesb$beta   = c(fit.bayesb$ETA[[1]]$b * fit.bayesb$ETA[[1]]$d)
  
  return (list(fit = fit.bayesb, t = t.bayesb[3],
               rsse = norm(y.test - X.test %*% fit.bayesb$beta - fit.bayesb$mu, '2')))
}


#'
#'
#'
#' varbvs
fit.varbvs = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run varbvs
  t.varbvs          = system.time(
    fit.varbvs       <- varbvs(X, Z = NULL, y, verbose = FALSE))
  
  return (list(fit = fit.varbvs, t = t.varbvs[3],
               rsse = norm(y.test - predict(fit.varbvs, X.test), '2')))
  
}

#'
#'
#'
#' varbvs2
fit.varbvs2 = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run varbvs
  t.varbvs          = system.time(
    fit.varbvs       <- varbvs(X, Z = NULL, y, verbose = FALSE, logodds = seq(-log10(p),1,length.out = 40)))
  
  return (list(fit = fit.varbvs, t = t.varbvs[3],
               rsse = norm(y.test - predict(fit.varbvs, X.test), '2')))
  
}

#'
#'
#'
#' susie
fit.susie = function(X, y, X.test, y.test, seed = 1, L = 20) {
  
  # set seed
  set.seed(seed)
  
  # run varbvs
  t.susie           = system.time(
    fit.susie        <- susie(X = X, Y = y, standardize = standardize, L = L))
  
  return (list(fit = fit.susie, t = t.susie[3],
               rsse = norm(y.test - predict(fit.susie, X.test), '2')))
  
}