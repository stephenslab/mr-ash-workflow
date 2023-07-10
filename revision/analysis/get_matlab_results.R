

n            = 500
p            = 200
s_range      = c(1,2,5,10,20,50,100,200)
trimmedlasso = matrix(0, 20, 8)

for (iter in 1:8) {
  s               = s_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_normal_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso[i,iter],"\n")
  }
}


trimmedlasso

tlasso = colMeans(trimmedlasso)



n            = 500
p            = 1000
s2_range      = c(1,3,10,30,100,300,1000)
trimmedlasso2 = matrix(0, 20, 7)

for (iter in 1:7) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_normal_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso2[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso2[i,iter],"\n")
  }
}

trimmedlasso2

tlasso2 = colMeans(trimmedlasso2)





n            = 500
p            = 1000
s2_range      = c(1,3,10,30,100,300,1000)
trimmedlasso3 = matrix(0, 20, 7)

for (iter in 1:7) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", design = "equicorrgauss", rho = 0.50, pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_rho050_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso3[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso3[i,iter],"\n")
  }
}

trimmedlasso3

tlasso3 = colMeans(trimmedlasso3)



aa = data.frame("t1" = tlasso4, row.names = c(5, 20 , 500))

n            = 500
p            = 1000
s2_range      = c(5, 20, 500, 1000)
trimmedlasso4 = matrix(0, 20, 3)
sparsity4 = matrix(0, 20, 3)

aa = data.frame(row.names = c(5, 20 , 500, 1000))

for (kk in 1:4) {
  for (iter in 1:4) {
    s               = s2_range[iter]
    for (i in 1:20) {
      data          = simulate_data(n, p, s = s, seed = i, signal = "t", df = 2^(kk-1), pve = 0.5)
      name = sprintf("result_fast_n%dp%ds%d_t%d_seed%d.mat", n,p, s, kk, i)
      name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
      b = readMat(name)$b
      sparsity4[i,iter] = sum(b != 0)
      trimmedlasso4[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
      cat(name,iter,i,s,trimmedlasso4[i,iter],"\n")
    }
  }
  aa[kk] = colMeans(trimmedlasso4)
}

for (iter in 1:4) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "t", df = 1, pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_t1_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    sparsity4[i,iter] = sum(b != 0)
    trimmedlasso4[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso4[i,iter],"\n")
  }
}

tlasso4 = colMeans(trimmedlasso4)
aa$const = tlasso4
aa

compute_oracle_ols <- function(data) {
  n = dim(data$X)[1]
  s = sum(data$beta != 0)
  if (s > n) {
    return (NA)
  }
  if (s == 0) {
    return (NA)
  }
  if (s == 1) {
    ind = which(data$beta != 0)
    #X = data$X[,ind]
    #X.test = data$X.test[,ind]
    #ols = lm(y ~ x, data.frame(y = data$y, x = X))
    #pred = predict(ols, data.frame(x = X.test))
    #err = norm(data$y.test - pred, '2') / data$sigma / sqrt(n)
    X = data$X[,c(ind, ind)]
    X.test = data$X.test[,c(ind, ind)]
    ols = glmnet(x = X, y = data$y, lambda = 0)
    pred = predict(ols, X.test)
    err = norm(data$y.test - pred, '2') / data$sigma / sqrt(n)
  } else {
    ind = (data$beta != 0)
    X = data$X[,ind]
    X.test = data$X.test[,ind]
    ols = glmnet(x = X, y = data$y, lambda = 0)
    pred = predict(ols, X.test)
    err = norm(data$y.test - pred, '2') / data$sigma / sqrt(n)
  }
  return (err)
}


################################################################
n            = 500
p            = 200
s_range      = c(1,2,5,10,20,50,100,200)
ols1          = matrix(0, 20, 8)
trimmedlasso1 = matrix(0, 20, 8)

for (iter in 1:8) {
  s               = s_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_normal_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso1[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    ols1[i,iter] = compute_oracle_ols(data)
    cat(name,iter,i,s,trimmedlasso1[i,iter],"\n")
  }
}
tlasso1 = colMeans(trimmedlasso1)
olsres1 = colMeans(ols1)
################################################################

################################################################
n            = 500
p            = 200
s_range      = c(1,2,5,10,20,50,100,200)
ols1          = matrix(0, 20, 8)
trimmedlasso1 = matrix(0, 20, 8)

for (iter in 1:8) {
  s               = s_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_normal_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso1[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    ols1[i,iter] = compute_oracle_ols(data)
    cat(name,iter,i,s,trimmedlasso1[i,iter],"\n")
  }
}
tlasso1 = colMeans(trimmedlasso1)
olsres1 = colMeans(ols1)
################################################################

################################################################
n            = 500
p            = 1000
s2_range      = c(1,3,10,30,100,300,1000)
ols2          = matrix(0, 20, 7)
trimmedlasso2 = matrix(0, 20, 7)

for (iter in 1:7) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_normal_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso2[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    ols2[i,iter] = compute_oracle_ols(data)
    cat(name,iter,i,s,trimmedlasso2[i,iter],"\n")
  }
}
tlasso2 = colMeans(trimmedlasso2)
olsres2 = colMeans(ols2)
################################################################

################################################################
n            = 500
p            = 1000
s2_range      = c(1,3,10,30,100,300,1000)
trimmedlasso3 = matrix(0, 20, 7)
ols3          = matrix(0, 20, 7)

for (iter in 1:7) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", design = "equicorrgauss", rho = 0.5, pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_rho050_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso3[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    ols3[i,iter] = compute_oracle_ols(data)
    cat(name,iter,i,s,trimmedlasso3[i,iter],"\n")
  }
}
tlasso3 = colMeans(trimmedlasso3)
olsres3 = colMeans(ols3)
################################################################

################################################################
n            = 500
p            = 1000
s2_range      = c(1,3,10,30,100,300,1000)
trimmedlasso5 = matrix(0, 20, 7)
ols5          = matrix(0, 20, 7)

for (iter in 1:7) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", design = "equicorrgauss", rho = 0.95, pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_rho095_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso5[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    ols5[i,iter] = compute_oracle_ols(data)
    cat(name,iter,i,s,trimmedlasso5[i,iter],"\n")
  }
}
tlasso5 = colMeans(trimmedlasso5)
olsres5 = colMeans(ols5)
################################################################

################################################################
n            = 500
p            = 10000
s2_range      = c(1,5,20,100,500,2000,10000)
trimmedlasso6 = matrix(0, 20, 7)
ols6          = matrix(0, 20, 7)

for (iter in 1:7) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = 0.5)
    #data          = readMat(sprintf("C:/Users/bions/Desktop/matdata/n%dp%ds%d_normal_seed%d.mat", n, p, s, i))
    name = sprintf("result_fast_n%dp%ds%d_normal_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols6[i,iter] = compute_oracle_ols(data)
    trimmedlasso6[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso6[i,iter],"\n")
  }
}
tlasso6 = colMeans(trimmedlasso6)
olsres6 = colMeans(ols6)
################################################################

################################################################
n            = 500
p            = 1000
s2_range      = c(5, 20, 500, 1000)
trimmedlasso7 = matrix(0, 20, 3)
sparsity7 = matrix(0, 20, 3)
ols5          = matrix(0, 20, 7)

for (iter in 1:4) {
  s               = s2_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "t", df = 1, pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_pve0.990000_seed%d.mat", n,p, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    sparsity7[i,iter] = sum(b != 0)
    trimmedlasso7[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso7[i,iter],"\n")
  }
}
################################################################



################################################################
n            = 500
p            = 1000
s = 5
trimmedlasso8 = matrix(0, 20, 8)
ols8          = matrix(0, 20, 8)
signalname = c("t", "t", "t", "t", "lap", "normal", "unif", "const")
for (kk in 1:8) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = signalname[kk], df = 2^(kk-1), pve = 0.5)
    if (kk <= 4){
      sname       = sprintf("t%s", kk)
    } else {
      sname       = signalname[kk]
    }
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols8[i,kk] = compute_oracle_ols(data)
    trimmedlasso8[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso8[i,kk],"\n")
  }
}
olsres8 = colMeans(ols8)
tlasso8 = colMeans(trimmedlasso8)
################################################################

################################################################
n = 500
p = 1000
s = 20
ols9          = matrix(0, 20, 8)
trimmedlasso9 = matrix(0, 20, 8)
signalname = c("t", "t", "t", "t", "lap", "normal", "unif", "const")
for (kk in 1:8) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = signalname[kk], df = 2^(kk-1), pve = 0.5)
    if (kk <= 4){
      sname       = sprintf("t%s", kk)
    } else {
      sname       = signalname[kk]
    }
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols9[i,kk] = compute_oracle_ols(data)
    trimmedlasso9[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso9[i,kk],"\n")
  }
}
olsres9 = colMeans(ols9)
tlasso9 = colMeans(trimmedlasso9)
################################################################

################################################################
n = 500
p = 1000
s = 500
ols10          = matrix(0, 20, 8)
trimmedlasso10 = matrix(0, 20, 8)
signalname = c("t", "t", "t", "t", "lap", "normal", "unif", "const")
for (kk in 1:8) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = signalname[kk], df = 2^(kk-1), pve = 0.5)
    if (kk <= 4){
      sname       = sprintf("t%s", kk)
    } else {
      sname       = signalname[kk]
    }
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols10[i,kk] = compute_oracle_ols(data)
    trimmedlasso10[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso10[i,kk],"\n")
  }
}
olsres10 = colMeans(ols10)
tlasso10 = colMeans(trimmedlasso10)
################################################################

################################################################
n = 500
p = 1000
s = 1000
ols11          = matrix(0, 20, 8)
trimmedlasso11 = matrix(0, 20, 8)
signalname = c("t", "t", "t", "t", "lap", "normal", "unif", "const")
for (kk in 1:8) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = signalname[kk], df = 2^(kk-1), pve = 0.5)
    if (kk <= 4){
      sname       = sprintf("t%s", kk)
    } else {
      sname       = signalname[kk]
    }
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols11[i,kk] = compute_oracle_ols(data)
    trimmedlasso11[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso11[i,kk],"\n")
  }
}
olsres11 = colMeans(ols11)
tlasso11 = colMeans(trimmedlasso11)
################################################################

################################################################
n = 500
p = 1000
s = 5
trimmedlasso12 = matrix(0, 20, 14)
ols12          = matrix(0, 20, 14)
pve_range    = c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
for (kk in 1:14) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = pve_range[kk])
    ppve          = pve_range[kk]
    sname = sprintf("pve%f", ppve)
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso12[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    ols12[i,kk] = compute_oracle_ols(data)
    cat(name,iter,i,s,trimmedlasso12[i,kk],"\n")
  }
}
tlasso12 = colMeans(trimmedlasso12)
olsres12 = colMeans(ols12)
################################################################

################################################################
s = 20
n = 500
p = 1000
trimmedlasso13 = matrix(0, 20, 14)
ols13          = matrix(0, 20, 14)
pve_range    = c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
for (kk in 1:14) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = pve_range[kk])
    ppve = pve_range[kk]
    sname = sprintf("pve%f", ppve)
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    trimmedlasso13[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    ols13[i,kk] = compute_oracle_ols(data)
    cat(name,iter,i,s,trimmedlasso13[i,kk],"\n")
  }
}
tlasso13 = colMeans(trimmedlasso13)
olsres13 = colMeans(ols13)
################################################################

################################################################
s = 500
n = 500
p = 1000
ols14          = matrix(0, 20, 14)
trimmedlasso14 = matrix(0, 20, 14)
pve_range    = c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
for (kk in 1:14) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = pve_range[kk])
    ppve = pve_range[kk]
    sname = sprintf("pve%f", ppve)
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols14[i,kk] = compute_oracle_ols(data)
    trimmedlasso14[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso14[i,kk],"\n")
  }
}
tlasso14 = colMeans(trimmedlasso14)
olsres14 = colMeans(ols14)
################################################################

################################################################
s = 1000
n = 500
p = 1000
ols15          = matrix(0, 20, 14)
trimmedlasso15 = matrix(0, 20, 14)
pve_range    = c(0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)
for (kk in 1:14) {
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = pve_range[kk])
    ppve = pve_range[kk]
    sname = sprintf("pve%f", ppve)
    name = sprintf("result_fast_n%dp%ds%d_%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols15[i,iter] = compute_oracle_ols(data)
    trimmedlasso15[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso15[i,kk],"\n")
  }
}
tlasso15 = colMeans(trimmedlasso15)
olsres15 = colMeans(ols15)
################################################################


################################################################
n = 287
s = 1000
trimmedlasso16 = matrix(0, 20, 5)
ols16 = matrix(0, 20, 5)
s2_range      = c(1, 5, 20, 100, 500)
filelist     = paste("data/", list.files("data", pattern = "*.RDS"), sep = "")
for (kk in 1:5) {
  s = s2_range[kk]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", design = "realgenotype", filepath = filelist[i])
    ppve = pve_range[kk]
    sname = sprintf("pve%f", ppve)
    name = sprintf("result_fast_n%ds%d_geno_seed%d.mat", n, s, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols16[i,kk] = compute_oracle_ols(data)
    trimmedlasso16[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso16[i,kk],"\n")
  }
}
tlasso16 = colMeans(trimmedlasso16)
olsres16 = colMeans(ols16)
################################################################

################################################################
n            = 500
p_range      = c(20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000)
s = 20
trimmedlasso17 = matrix(0, 20, 10)
tritime17 = matrix(0, 20, 10)
ols17 = matrix(0, 20, 10)

for (iter in 1:10) {
  p               = p_range[iter]
  for (i in 1:20) {
    data          = simulate_data(n, p, s = s, seed = i, signal = "normal", pve = pve_range[kk])
    name = sprintf("result_fast_n%dp%ds20_largenfixedp_seed%d.mat", n, p, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols17[i,iter] = compute_oracle_ols(data)
    trimmedlasso17[i,iter] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    name = sprintf("time_fast_n%dp%ds20_largenfixedp_seed%d.mat", n, p, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    tritime17[i, iter] = readMat(name)$t
    cat(name,iter,i,s,trimmedlasso17[i,iter],"\n")
  }
}
tlasso17 = colMeans(trimmedlasso17)
olsres17 = colMeans(ols17)
ttime17 = colMeans(tritime17)
tmed17 = apply(tritime17, 2, median)
################################################################


################################################################
n = 500
p = 1000
s = 20
ols18          = matrix(0, 20, 7)
trimmedlasso18 = matrix(0, 20, 7)
signalname = c("t1", "t2", "t4", "t8", "lap", "normal", "unif", "const")
for (kk in 1:7) {
  for (i in 1:20) {
    sname       = signalname[kk]
    data          = simulate_data(n, p, s = s, seed = i, signal = 'normal', noise = sname, df = 2^(kk-1), pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_noise%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols18[i,kk] = compute_oracle_ols(data)
    trimmedlasso18[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso18[i,kk],"\n")
  }
}
olsres18 = colMeans(ols18)
tlasso18 = colMeans(trimmedlasso18)
################################################################

################################################################
n = 500
p = 1000
s = 500
ols19          = matrix(0, 20, 7)
trimmedlasso19 = matrix(0, 20, 7)
signalname = c("t1", "t2", "t4", "t8", "lap", "normal", "unif", "const")
for (kk in 1:7) {
  for (i in 1:20) {
    sname       = signalname[kk]
    data          = simulate_data(n, p, s = s, seed = i, signal = 'normal', noise = sname, df = 2^(kk-1), pve = 0.5)
    name = sprintf("result_fast_n%dp%ds%d_noise%s_seed%d.mat", n,p, s, sname, i)
    name = paste("C:/Users/bions/Documents/MATLAB/", name, sep = "")
    b = readMat(name)$b
    ols19[i,kk] = compute_oracle_ols(data)
    trimmedlasso19[i,kk] = norm(data$y.test - data$X.test %*% as.vector(b), '2') / data$sigma / sqrt(n)
    cat(name,iter,i,s,trimmedlasso19[i,kk],"\n")
  }
}
olsres19 = colMeans(ols19)
tlasso19 = colMeans(trimmedlasso19)
################################################################