# TO DO: Explain here what this script does, and how to use it.
library(glmnet)

# SIMULATE DATA
# -------------
set.seed(1)
n <- 500
p <- 2000
X <- matrix(rnorm(n*p),n,p)

# Generate regression coefficients so that only 20 of the predictors
# affect the response.
i    <- sample(p,20)
b    <- rep(0,p)
b[i] <- rnorm(na)

# Generate the responses.
y <- drop(X %*% b + rnorm(n))

# FIT ELASTIC NET MODEL
# ---------------------
# TO DO.

# FIT BAYESR MODEL
# ----------------
# TO DO.
