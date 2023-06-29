# A small experiment to compare performance of least-squares (OLS)
# estimates vs. mr.ash in small, simulated data sets. This experiment
# supports our response to Comment 5 from Reviewer 1.
library(glmnet)
library(mr.ash.alpha)
library(ggplot2)
library(cowplot)
n <- 200
p <- c(2,4,8,16,32,64)
r <- 0.5
ns <- 20
np <- length(p)
set.seed(1)

# Set up the data structures for storing the results.
rmse_glmnet <- matrix(0,np,ns)
rmse_mrash  <- matrix(0,np,ns)

# Repeat for each setting of p, the number of predictors.
for (i in 1:np) {

  # Repeat for each simulation.
  for (j in 1:ns) {
    cat(sprintf("p = %d, simulation %d\n",p[i],j))
    
    # Simulate a training and test data set.
    X  <- matrix(rnorm(2*n*p[i]),2*n,p[i])
    X  <- scale(X,center = TRUE,scale = TRUE)
    b  <- rnorm(p[i])
    sb <- r/(1-r)/var(drop(X %*% b))
    b  <- sqrt(sb) * b
    y  <- drop(X %*% b + rnorm(n))
    y  <- y - mean(y)
    rows  <- seq(n+1,2*n)
    Xtest <- X[rows,]
    ytest <- y[rows]
    rows  <- 1:n
    X <- X[rows,]
    y <- y[rows]

    # Compute least-squares estimates of the coefficients.
    fit_glmnet <- glmnet(X,y,family = "gaussian",lambda = 0,
                         standardize = TRUE,intercept = TRUE)

    # Compute the glmnet test set predictions.
    ypred <- drop(predict(fit_glmnet,Xtest))

    # Compute the RMSE.
    rmse_glmnet[i,j] <- sqrt(mean((ytest - ypred)^2))
    
    # Fit a mr.ash model.
    capture.output(suppressWarnings(
      fit_mrash <- mr.ash(X,y,standardize = TRUE,intercept = TRUE)))

    # Compute the mr.ash test set predictions.
    ypred <- drop(predict(fit_mrash,Xtest))
    
    # Compute the RMSE.
    rmse_mrash[i,j] <- sqrt(mean((ytest - ypred)^2))
  }
}

rmse_null <- 1/sqrt(0.5)

# Summarize the results in a plot.
pdat <- data.frame(p   = as.vector(matrix(p,np,ns,byrow = FALSE)),
                   sim = as.vector(matrix(1:ns,np,ns,byrow = TRUE)),
                   y   = as.vector((rmse_mrash - rmse_glmnet)/rmse_null),
                   stringsAsFactors = FALSE)
pdat <- transform(pdat,p = factor(p))
p1 <- ggplot(pdat,aes(x = p,y = y)) +
  geom_boxplot(width = 0.5) +
  labs(x = "number of predictors",y = "scaled RMSE mr.ash - OLS") +
  theme_cowplot(font_size = 12)
ggsave("mrash_vs_ls.eps",p1,height = 2.75,width = 3)
