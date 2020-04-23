# TO DO: Explain here what this script does, and how to use it.
library(glmnet)

# SIMULATE DATA
# -------------
set.seed(1)
n <- 500
p <- 2000
X <- matrix(rnorm(n*p),n,p)
X <- scale(X,center = TRUE,scale = FALSE)

# Generate regression coefficients so that only 5 of the predictors
# affect the response.
i    <- sample(p,5)
b    <- rep(0,p)
b[i] <- rnorm(5)

# Generate the responses.
y <- drop(X %*% b + 2*rnorm(n))

# FIT ELASTIC NET MODEL
# ---------------------
# Fit the elastic net model, and use it to estimate the responses in
# the training data.
fit.glmnet <- cv.glmnet(X,y,alpha = 0.95)
yest <- drop(predict(fit.glmnet,X,s = fit.glmnet$lambda.min))
plot(y,yest,pch = 20,col = "darkblue")
abline(a = 0,b = 1,lty = "dotted",col = "magenta")

# FIT BAYESR MODEL
# ----------------
# TO DO.
