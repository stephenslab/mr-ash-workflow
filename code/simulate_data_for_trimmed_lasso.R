# module load  R/3.5.1
library(R.matlab)
source("sim_wrapper.R")
n <- 500
p <- 1e4
s_range <- c(1,5,20,100,500,2000,10000)
for (iter in 1:7) {
  s <- s_range[iter]
  for (i in 1:20) {
    cat("s = ",s,", i = ",i,"\n",sep="")
    dat <- simulate_data(n,p,s = s,seed = i,signal = "normal",pve = 0.5)
    fname <- sprintf("n%dp%ds%d_normal_seed%d.mat",n,p,s,i)
    writeMat(fname,X = dat$X,y = dat$y)
  }
}
