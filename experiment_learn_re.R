# Which method is better at learning random effect?
rm(list = ls())
suppressMessages(require(mvtnorm))
# Simulation
p = 10
n = 40
raneff = rnorm(p)
Sfct = matrix(0.5*rnorm(2*p),p,2)
Sigma = Sfct%*%t(Sfct)+diag(0.5,p)
avs = t(rmvnorm(n,-3+raneff,Sigma))
y = matrix(rpois(n*p,exp(c(avs))),p,n)

# test function
test <- function(y,mu,Sigma,nm,method){
  suppressMessages( require(Rcpp) )
  source("~/Dropbox/PhD Research/Hierarchical Random Effect Model_Statistical Ecology/gmm_ars/ars_fcn.R")
  sourceCpp("~/Dropbox/PhD Research/Hierarchical Random Effect Model_Statistical Ecology/Eric/updateBeta.cpp")
  sd_re = 100
  Sigma_inv = solve(Sigma)
  Rinfo <- condNorm(Sigma)
  p = nrow(y)
  n = ncol(y)
  vre = NULL
  re = rep(0,p)
  av = matrix(0,p,n)
  for (iter in 1:nm){
    if (iter%%floor(nm/10)==0){
      cat(iter,"\n")
    }
    for (i in 1:n){
      av[,i] = cars_c(y[,i], av[,i], mu+re, Rinfo)
    }
    tmp = mMmu_M(mu, av, FALSE)
    tmpS = solve(n*Sigma_inv+diag(1/sd_re,p))
    re = rmvnorm(1,tmpS%*%apply(Sigma_inv%*%tmp,1,sum),tmpS)
    vre = cbind(vre,c(re))
  }
  return(vre)
}

# experiment
nm = 10000
res = test(y,rep(-3,p),Sigma,nm,"pm")

