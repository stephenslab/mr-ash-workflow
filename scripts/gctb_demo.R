# TO DO: Explain here what this script does, and how to use it.
library(glmnet)
source("../code/gctb.R")

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

# FIT BayesR MODEL
# ----------------
# Create the "fam" file.
dat <- data.frame(fid   = paste0("s",1:n),
                  iid   = paste0("s",1:n),
                  pid   = 0,
                  mid   = 0,
                  sex   = 0,
                  pheno = -9)
write.table(dat,"sim.fam",sep = " ",row.names = FALSE,col.names = FALSE,
            quote = FALSE)

# Create the "map" file.
dat <- data.frame(chr = 1,
                  id  = paste0("rs",1:p),
                  cM  = 0,
                  bp  = 1000*(1:p))
write.table(dat,"sim.bim",sep = " ",row.names = FALSE,col.names = FALSE,
            quote = FALSE)

# Create the "ped" file.


# Create the PLINK bed/bim/fam fileset.


#
#   gctb --bfile 1000G_eur_chr22 --pheno test.phen \
#     --bayes R --seed 1 --gamma 0,0.01,0.1,1 --pi 0.25,0.25,0.25,0.25
#     --chain-length 1500 --burn-in 500 --out test 
#
