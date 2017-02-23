# Sample the auxillary variables in Log-Gaussian Mixture of Poisson

# Load dependencies
suppressMessages(require(mvtnorm))
suppressMessages( require(Rcpp) )
sourceCpp("~/Dropbox/PhD Research/Hierarchical Random Effect Model_Statistical Ecology/Eric/update_mu.cpp")

# Helper Functions
# Method 1. Independent Metropolis Hastings
imh <- function(y,mu,Sigma,x){
  p = length(y)
  nx = c(rmvnorm(1,mu,Sigma))
  lar = sum((nx-x)*y)-sum(exp(nx))+sum(exp(x))
  if (log(runif(1))<lar){
    x = nx
  }
  return(list(newx = x,
         lar = lar))
}

# Method 2. Pure Metropolis
pm <- function(y,mu,Sigma,x){
  p = length(y)
  nx = c(rmvnorm(1,x,Sigma))
  lar = sum((nx-x)*y)-sum(exp(nx))+sum(exp(x))+dmvnorm(nx,mu,Sigma,log=TRUE)
  lar = lar - dmvnorm(x,mu,Sigma,log=TRUE)
  if (log(runif(1))<lar){
    x = nx
  }
  return(list(newx = x,
              lar = lar))
}

# Main function
sampleLGP <- function(y,mu,Sigma,nm,method="conARS"){
  if (method=="conARS"){
    Rinfo <- condNorm(Sigma)
  }
  p = nrow(y)
  x = rep(0,p)
  sps = NULL
  lars = NULL
  for (iter in 1:nm){
    if (method=="imh"){
      tmp = imh(y,mu,Sigma,x)
      sps = cbind(sps,tmp$newx)
      lars = c(lars,tmp$lar)
    }else if (method=="pm"){
      tmp = pm(y,mu,Sigma,x)
      sps = cbind(sps,tmp$newx)
      lars = c(lars,tmp$lar)
    }else{
      x = cars_c(y, x, mu, Rinfo)
      sps = cbind(sps,x)
    }
  }
  return(list(samples = sps,
              lars = lars))
}
