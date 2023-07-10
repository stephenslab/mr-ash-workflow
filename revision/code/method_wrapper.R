#'
if (exists("standardize")) {
  if (is.logical(standardize)) {
    
  } else {
    warning("standardize must be logical: set to FALSE")
    standardize     = FALSE
  }
} else {
  standardize       = FALSE
}

#'
#'
#'
get_phi <- function(fit) {
  # compute residual
  r            = fit$data$y - fit$data$X %*% fit$beta
  
  # compute bw and S2inv
  bw           = as.vector((t(fit$data$X) %*% r) + fit$data$w * fit$beta)
  S2inv        = 1 / outer(fit$data$w, 1/fit$data$sa2, '+');
  
  # compute mu, phi
  mu           = bw * S2inv;
  phi          = -log(1 + outer(fit$data$w, fit$data$sa2))/2 + mu * (bw / 2 / fit$sigma2);
  phi          = c(fit$pi) * t(exp(phi - apply(phi,1,max)));
  phi          = t(phi) / colSums(phi);
  return (list(phi = phi, mu = mu, r = r))
}

#'
#'
#'
#' MR.ASH
fit.mr.ash = function(X, y, X.test, y.test, seed = 1, sa2 = NULL, ...) {
  
  # set seed
  set.seed(seed)
  
  # run mr.ash
  t.mr.ash           = system.time(
    fit.mr.ash        <- mr.ash(X = X, y = y, sa2 = sa2,
                                max.iter = 2000,
                                standardize = standardize,
                                tol = list(epstol = 1e-12, convtol = 1e-8), ...))
  beta               = fit.mr.ash$beta
  pip                = 1 - get_phi(fit.mr.ash)$phi[,1]
  
  return (list(fit = fit.mr.ash, t = t.mr.ash[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict(fit.mr.ash, X.test), '2')))
}

#'
#'
#'
#' MR.ASH different sigma2 update
fit.mr.ash.sigma2 = function(X, y, X.test, y.test, beta.init = NULL, update.order = NULL,
                             seed = 1, sa2 = NULL, ...) {
  
  # set seed
  set.seed(seed)
  
  # run mr.ash
  t.mr.ash           = system.time(
    fit.mr.ash        <- mr.ash(X = X, y = y, sa2 = sa2, method = "sigma",
                                max.iter = 2000, beta.init = beta.init,
                                update.order = update.order,
                                standardize = standardize,
                                tol = list(epstol = 1e-12, convtol = 1e-8), ...))
  beta               = fit.mr.ash$beta
  pip                = 1 - get_phi(fit.mr.ash)$phi[,1]
  
  return (list(fit = fit.mr.ash, t = t.mr.ash[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict(fit.mr.ash, X.test), '2')))
}

#'
#'
#'
#' MR.ASH
fit.mr.ash2 = function(X, y, X.test, y.test, seed = 1,
                       sa2 = NULL,
                       update.pi = TRUE, pi = NULL,
                       beta.init = NULL, sigma2 = NULL,
                       update.order = NULL, ...) {
  
  
  # set seed
  set.seed(seed)
  
  # run mr.ash order
  t.mr.ash           = system.time(
    fit.mr.ash        <- mr.ash(X = X, y = y, sa2 = sa2, update.order = update.order,
                                max.iter = 2000, min.iter = 200,
                                beta.init = beta.init, update.pi = update.pi, pi = pi,
                                standardize = standardize, sigma2 = sigma2,
                                tol = list(epstol = 1e-12, convtol = 1e-8), ...))
  beta               = fit.mr.ash$beta
  pip                = 1 - get_phi(fit.mr.ash)$phi[,1]
  
  return (list(fit = fit.mr.ash, t = t.mr.ash[3], beta = beta, pip = pip,
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
  beta              = as.vector(coef(fit.lasso, s = fit.lasso$lambda.1se))[-1]
  beta2             = as.vector(coef(fit.lasso, s = fit.lasso$lambda.min))[-1]
  
  return (list(fit = fit.lasso, t = t.lasso[3], t2 = t.lasso2[3], beta = beta, beta2 = beta2,
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
  beta              = as.vector(coef(fit.ridge, s = fit.ridge$lambda.1se))[-1]
  beta2             = as.vector(coef(fit.ridge, s = fit.ridge$lambda.min))[-1]
  
  return (list(fit = fit.ridge, t = t.ridge[3], beta = beta, beta2 = beta2,
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
  
  beta              = as.vector(coef(fit.enet[[alpha.ind]], s = fit.enet[[alpha.ind]]$lambda.1se))[-1]
  beta2             = as.vector(coef(fit.enet[[alpha.ind]], s = fit.enet[[alpha.ind]]$lambda.min))[-1]
  
  return (list(fit = fit.enet[[alpha.ind]], t = sum(time), beta = beta, beta2 = beta2,
               rsse = rsse[alpha.ind], rsse2 = rsse2[alpha.ind], alpha = (alpha.ind - 1) * 0.1))
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
  beta              = as.vector(coef(fit.scad))[-1]
  
  return (list(fit = fit.scad, t = t.scad[3], beta = beta, rsse = norm(y.test - predict(fit.scad, X.test), '2')))
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
  beta              = as.vector(coef(fit.mcp))[-1]
  
  return (list(fit = fit.mcp, t = t.mcp[3], beta = beta, rsse = norm(y.test - predict(fit.mcp, X.test), '2')))
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
  beta              = as.vector(coef(fit.scad))[-1]
  
  return (list(fit = fit.scad, t = t.scad[3], beta = beta, rsse = norm(y.test - predict(fit.scad, X.test), '2'),
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
  t.mcp             = system.time(
    fit.mcp        <- cv.ncvreg(X, y, penalty = "MCP", nfolds = 10))
  t.mcp2            = system.time(ncvreg(X, y, penalty = "MCP"))
  beta              = as.vector(coef(fit.mcp))[-1]
  
  return (list(fit = fit.mcp, t = t.mcp[3], beta = beta, rsse = norm(y.test - predict(fit.mcp, X.test), '2'),
               t2 = t.mcp2[3]))
}


#'
#'
#'
#' L0Learn
fit.l0l2learn = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run l0learn
  t.l0learn         = system.time(
    fit.l0learn      <- L0Learn.cvfit(X, y, nFolds = 10, penalty = "L0L2"))

  cvres             = print(fit.l0learn)
  cvres$error       = unlist(fit.l0learn$cvMeans)
  cvmin             = cvres[which.min(cvres$error),]
  lambda.min        = cvmin$lambda
  gamma.min         = cvmin$gamma
  beta              = as.vector(coef(fit.l0learn, lambda = lambda.min, gamma = gamma.min))[-1]
  
  return (list(fit = fit.l0learn, t = t.l0learn[3], lambda = lambda.min, beta = beta,
               rsse = norm(y.test - as.vector(predict(fit.l0learn, X.test, lambda = lambda.min, gamma = gamma.min)[,1]), '2')))
}

#'
#'
#'
#' L0Learn
fit.l0l1learn = function(X, y, X.test, y.test, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # run l0learn
  t.l0learn         = system.time(
    fit.l0learn      <- L0Learn.cvfit(X, y, nFolds = 10, penalty = "L0L1"))
  
  cvres             = print(fit.l0learn)
  cvres$error       = unlist(fit.l0learn$cvMeans)
  cvmin             = cvres[which.min(cvres$error),]
  lambda.min        = cvmin$lambda
  gamma.min         = cvmin$gamma
  beta              = as.vector(coef(fit.l0learn, lambda = lambda.min, gamma = gamma.min))[-1]

  return (list(fit = fit.l0learn, t = t.l0learn[3], lambda = lambda.min, gamma = gamma.min, beta = beta,
               rsse = norm(y.test - as.vector(predict(fit.l0learn, X.test, lambda = lambda.min, gamma = gamma.min)[,1]), '2')))
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
  beta              = as.vector(coef(fit.l0learn, lambda = lambda.min))[-1]
  
  return (list(fit = fit.l0learn, t = t.l0learn[3], lambda = lambda.min, beta = beta,
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
  
  return (list(fit = fit.blasso, t = t.blasso[3], beta = fit.blasso$beta,
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
  pip               = (fit.bayesb$ETA[[1]]$d > 0.5)
  
  return (list(fit = fit.bayesb, t = t.bayesb[3], beta = fit.bayesb$beta, pip = pip,
               rsse = norm(y.test - X.test %*% fit.bayesb$beta - fit.bayesb$mu, '2')))
}

#'
#'
#'
#' BayesC
fit.bayesc = function(X, y, X.test, y.test, seed = 1, nIter = NULL, burnIn = NULL) {
  
  # set seed
  set.seed(seed)
  
  # default setting
  if (is.null(nIter)) nIter = 1500
  if (is.null(burnIn)) burnIn = 500
  
  # run bayesb
  t.bayesc          = system.time(
    fit.bayesc       <- BGLR(y, ETA = list(list(X = X, model="BayesC", standardize = standardize)),
                             verbose = FALSE, nIter = nIter, burnIn = burnIn))
  fit.bayesc$beta   = c(fit.bayesc$ETA[[1]]$b * fit.bayesc$ETA[[1]]$d)
  
  return (list(fit = fit.bayesc, t = t.bayesc[3], beta = fit.bayesc$beta,
               rsse = norm(y.test - X.test %*% fit.bayesc$beta - fit.bayesc$mu, '2')))
}

#'
#'
#'
#' gibbs sampling
fit.mcmc = function(X, y, X.test, y.test, seed = 1, nIter = NULL, burnIn = NULL, beta.init = data$beta + 0) {
  
  # set seed
  set.seed(seed)
  
  # default setting
  if (is.null(nIter)) nIter = 1500
  if (is.null(burnIn)) burnIn = 500
  
  # run bayesb
  t.mcmc           = system.time(
    fit.mcmc <- gibbs.sampling(data$X, data$y, sa2 = c(0, 1 / s), burn.in = burnIn, max.iter = nIter,
                               pi = c(1 - s/p, s/p), beta.init = beta.init, sigma2 = data$sigma^2))
  fit.mcmc$data$sa2 = c(0, 1 / s)
  fit.mcmc$pi       = pi = c(1 - s/p, s/p)
  pip              = 1 - get_phi(fit.mcmc)$phi[,1]
  
  return (list(fit = fit.mcmc, t = t.mcmc[3], beta = fit.mcmc$beta, pip = pip, 
               rsse = norm(y.test - X.test %*% fit.mcmc$beta - fit.mcmc$mu, '2')))
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
  beta              = c(rowSums(fit.varbvs$alpha * fit.varbvs$mu))
  pip               = fit.varbvs$pip
  
  return (list(fit = fit.varbvs, t = t.varbvs[3], beta = beta, pip = pip,
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
  beta              = c(rowSums(fit.varbvs$alpha * fit.varbvs$mu))
  pip               = fit.varbvs$pip
  
  return (list(fit = fit.varbvs, t = t.varbvs[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict(fit.varbvs, X.test), '2')))
  
}

#'
#'
#'
#' susie
fit.susie = function(X, y, X.test, y.test, seed = 1, L = 20) {
  
  # set seed
  set.seed(seed)
  
  # run susie
  t.susie           = system.time(
    fit.susie        <- susie(X = X, Y = y, standardize = standardize, L = L))
  beta              = coef(fit.susie)[-1]
  pip               = fit.susie$pip
  
  return (list(fit = fit.susie, t = t.susie[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict(fit.susie, X.test), '2')))
  
}

#'
#'
#'
#' emvs
fit.emvs = function(X, y, X.test, y.test, seed = 1, independent = FALSE, epsilon =  1e-10,
                    v0 = exp(seq(-20, -1, length.out = 20)), v1 = 1000) {
  
  # set seed
  set.seed(seed)
  
  # run emvs
  t.emvs            = system.time(temp <- capture.output(
    fit.emvs        <- EMVS(Y = y, X = X, standardize = standardize, independent = independent,
                            v0 = v0, v1 = v1, epsilon = epsilon)))
  
  beta              = coef.emvs(fit.emvs)
  pip               = pip.emvs(fit.emvs)
  
  return (list(fit = fit.emvs, t = t.emvs[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict.emvs(fit.emvs, X.test), '2')))
  
}
#'
#'
#'
#' coef.emvs
coef.emvs = function(fit) {
  best_model        = round(median(which(fit$log_g_function == EMVSbest(fit)$log_g_function)))
  beta              = fit$betas[best_model,]
  return (beta)
}
pip.emvs = function(fit) {
  best_model        = round(median(which(fit$log_g_function == EMVSbest(fit)$log_g_function)))
  pip               = fit$prob_inclusion[best_model,]
  return (pip)
}
predict.emvs = function(fit, X) {
  best_model        = round(median(which(fit$log_g_function == EMVSbest(fit)$log_g_function)))
  beta              = fit$beta[best_model,]
  return (X %*% beta)
}

#'
#'
#'
#' emvscv
fit.emvscv = function(X, y, X.test, y.test, seed = 1, independent = FALSE, epsilon =  1e-10,
                      v0 = exp(seq(-20, -1, length.out = 10)), v1 = 1000) {
  
  # set seed
  set.seed(seed)
  
  # cv prepare
  nv                      = 10
  n                       = nrow(X)
  m                       = as.integer(n / 5)
  cv_shuffle              = sample(n,n)
  cv_partition            = list(cv_shuffle[1:m],
                                 cv_shuffle[(m+1):(2*m)], 
                                 cv_shuffle[(2*m+1):(3*m)], 
                                 cv_shuffle[(3*m+1):(4*m)], 
                                 cv_shuffle[(4*m+1):n])
  cv_error                = rep(0, nv)
  cv_time                 = 0
  
  # run cv
  for (cv_ind in 1:5) {
    t.cvem               = system.time({temp <- capture.output(
      fit.cvem          <- EMVS(Y = y[-cv_partition[[cv_ind]]], X = X[-cv_partition[[cv_ind]], ],
                                standardize = standardize, independent = independent,
                                v0 = v0, v1 = v1, epsilon = epsilon))
    resid_cv_v0         <- X[cv_partition[[cv_ind]], ] %*% t(fit.cvem$betas) -  y[cv_partition[[cv_ind]]]
    rmse_cv_v0          <- apply(resid_cv_v0, 2, norm, '2')
    })
    cv_error            = cv_error + rmse_cv_v0^2
    cv_time             = cv_time + t.cvem[3]
  }
  
  # evaluate cv
  ind_best              = which(cv_error == min(cv_error))[1]
  v0_best               = v0[ind_best]
  
  t.emvscv              = system.time(temp <- capture.output(
    fit.emvscv         <- EMVS(Y = y, X = X, standardize = standardize, independent = independent,
                              v0 = v0_best, v1 = v1, epsilon = epsilon)))
  
  cv_time               = cv_time + t.emvscv[3]
  beta                  = fit.emvscv$beta[1,]
  pip                   = fit.emvscv$prob_inclusion[1,]
  
  return (list(fit = fit.emvs, t = cv_time, beta = beta, pip = pip,
               rsse = norm(y.test - X.test %*% beta, '2')))
  
}

#'
#'
#'
#' sslasso
fit.sslasso = function(X, y, X.test, y.test, seed = 1, penalty = "adaptive", eps =  1e-10,
                       lambda0 = NULL, lambda1 = NULL) {
  
  # set seed
  set.seed(seed)
  
  # run sslasso
  t.sslasso         = system.time(temp <- capture.output(
    fit.sslasso     <- SSLASSO(X = X, y = y, penalty = penalty, 
                               eps = eps)))
  
  beta              = coef.sslasso(fit.sslasso)
  pip               = NULL
  
  return (list(fit = fit.sslasso, t = t.sslasso[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict.sslasso(fit.sslasso, X.test), '2')))
  
}

#'
#'
#'
#' coef.sslasso
coef.sslasso = function(fit) {
  beta              = fit$beta[,length(fit$lambda0)]
  return (beta)
}
predict.sslasso = function(fit, X) {
  beta              = fit$beta[,length(fit$lambda0)]
  return (X %*% beta)
}

#'
#'
#'
#' sslasso with cv
fit.sslassocv = function(X, y, X.test, y.test, seed = 1, penalty = "adaptive", eps =  1e-10,
                         lambda0 = NULL, lambda1 = NULL) {
  
  # set seed
  set.seed(seed)
  
  # cv prepare
  nlambda                 = 100
  n                       = nrow(X)
  m                       = as.integer(n / 5)
  cv_shuffle              = sample(n,n)
  cv_partition            = list(cv_shuffle[1:m],
                                 cv_shuffle[(m+1):(2*m)], 
                                 cv_shuffle[(2*m+1):(3*m)], 
                                 cv_shuffle[(3*m+1):(4*m)], 
                                 cv_shuffle[(4*m+1):n])
  lambda0                 = seq(1, n, length.out = nlambda)
  cv_error                = rep(0, nlambda)
  cv_time                 = 0
  
  # run cv
  for (cv_ind in 1:5) {
      t.cvss              = system.time({temp <- capture.output(
        fit.cvss          <- SSLASSO(X = X[-cv_partition[[cv_ind]], ], y = y[-cv_partition[[cv_ind]]],
                                  penalty = penalty, lambda0 = lambda0, eps = eps))
        resid_cv_lambda0  <- X[cv_partition[[cv_ind]], ] %*% fit.cvss$beta -  y[cv_partition[[cv_ind]]]
        rmse_cv_lambda0   <- apply(resid_cv_lambda0, 2, norm, '2')
      })
      cv_error            = cv_error + rmse_cv_lambda0^2
      cv_time             = cv_time + t.cvss[3]
  }
  ind_best                = which(cv_error == min(cv_error))[1]
  lambda0_best            = lambda0[ind_best]
  
  # run sslasso
  t.sslasso               = system.time(temp <- capture.output(
    fit.sslasso           <- SSLASSO(X = X, y = y, penalty = penalty, lambda0 = lambda0_best,
                                     lambda1 = 1, eps = eps)))
  
  cv_time                 = cv_time + t.sslasso[3]
  beta                    = fit.sslasso$beta[,1]
  pip                     = NULL
  
  return (list(fit = fit.sslasso, t = cv_time, beta = beta, pip = pip,
               rsse = norm(y.test - X.test %*% beta, '2')))
  
}

#'
#'
#'
#' bayeslm - horseshoe
fit.bayeslm = function(X, y, X.test, y.test, seed = 1, prior = "horseshoe", burnin = 0, N = 100,
                       icept = TRUE, singular = TRUE) {
  
  # set seed
  set.seed(seed)
  
  # run sslasso
  t.bayeslm         = system.time(temp <- capture.output(
    fit.bayeslm     <- bayeslm(Y = y, X = X, prior = prior, N = N, burnin = burnin, singular = singular,
                               standardize = standardize, icept = icept)))
  beta              = coef.bayeslm(fit.bayeslm)
  pip               = NULL
  
  return (list(fit = fit.bayeslm, t = t.bayeslm[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict.bayeslm(fit.bayeslm, X.test), '2')))
  
}
#'
#'
#'
#' coef.bayeslm
coef.bayeslm = function(fit) {
  beta              = as.vector(colMeans(fit$beta))
  return (beta[-1])
}
predict.bayeslm = function(fit, X) {
  beta              = as.vector(colMeans(fit$beta))
  return (X %*% beta[-1] + beta[1])
}

#'
#'
#'
#' horseshoe
fit.horseshoe = function(X, y, X.test, y.test, seed = 1, method.tau = "halfCauchy", burn = 0, nmc = 100,
                         method.sigma = "Jeffreys") {
  
  # set seed
  set.seed(seed)
  
  # run sslasso
  t.horseshoe          = system.time(temp <- capture.output(
    fit.horseshoe     <- horseshoe(y = y, X = X, burn = burn, nmc = nmc,
                                   method.tau = method.tau, method.sigma = method.sigma)))
  
  beta                = coef.horseshoe(fit.horseshoe)
  pip                 = NULL
  
  return (list(fit = fit.horseshoe, t = t.horseshoe[3], beta = beta, pip = pip,
               rsse = norm(y.test - predict.horseshoe(fit.horseshoe, X.test), '2')))
  
}
#'
#'
#'
#' coef.horseshoe
coef.horseshoe = function(fit) {
  beta              = as.vector(fit$BetaHat)
  return (beta)
}
predict.horseshoe = function(fit, X) {
  beta              = as.vector(fit$BetaHat)
  return (X %*% beta)
}