# Metropolis Hastings witin Gibbs + Moment Matching
# Use condtional ARS to sample the auxillary variables
# Empirically update mu and Sigma using moment matching

rm(list = ls())

# dependencies
cat("loading dependencies...")
suppressMessages(require(parallel))
suppressMessages(require(mvtnorm))
suppressMessages(require(MCMCpack))
suppressMessages( require(Rcpp) )
source("~/Dropbox/MyGitbub/poisson-log-gaussian-solver/ars_fcn.R")
sourceCpp("~/Dropbox/PhD Research/Hierarchical Random Effect Model_Statistical Ecology/HMCEM/monte_carlo_int.cpp")
cat("finished\n")

# simulate data
cat("simulating the data...")
n = 5000
p = 3
mu_true = -3+0.5*rnorm(p)
Sigma_true = diag(1+0.5*rnorm(p)^2)
nABC = 200
hashABC = sample(1:nABC,n,replace = TRUE)
SigAbc = matrix(c(1,-0.8,0,-0.8,1,0,0,0,1),3,3)
# tmp = cbind(c(1,1,1),c(0.2,0.8,0.2))
# SigAbc = tmp%*%t(tmp)+SigAbc
MuAbc = t(rmvnorm(nABC,rep(0,p),SigAbc))
y = matrix(0,p,n)
for (i in 1:n){
  b = rmvnorm(1,mu_true+MuAbc[,hashABC[i]],Sigma_true)
  y[,i] = rpois(p,exp(b))
}
cat("finished\n")

nb = 6000
ns = 4000
nmm = 100
# 
# 1. preparation
# 1.1. preparation for parallel computing
# put data into a list for parallel proccessing
y_listABC = list()
for (j in 1:nABC){
  y_listABC[[j]] = y[,hashABC==j]
}
ncores = detectCores()-1 # number of cores to use
# 1.2. initial estimate of mu and Sigma
p = nrow(y)
n = ncol(y)
muABC = matrix(0,p,nABC)
res = mmatch(y, muABC, hashABC, nABC)
mu = res$mu
Sigma = res$Sigma
# mu = mu_true
# Sigma = Sigma_true
Sigma_inv = solve(Sigma)

# 1.3. preparation for MH
av = matrix(0,p,n) # auxillary variables
SigmaABC = diag(1,p)
SigmaABC_inv = solve(SigmaABC)
mu_mm = rep(0,p)
Sigma_mm = matrix(0,p,p)

# 1.4. store quantities of interest
vSigma = NULL
vmu = NULL
vS = NULL



# iteration
for (iter in 1:(nb+ns)){
  if (iter%%((nb+ns)/100)==0){
    cat("iteration",iter,"\n")
  }
  # update auxillary variables av
  Rinfo = condNorm(Sigma)
  tmp = mclapply(1:nABC,FUN = function(x){
    res = mcars_c(y_listABC[[x]], av[,hashABC==x], mu+muABC[,x], Rinfo)
    dmx = mMmu_M(mu, res, FALSE)
    tmpS = solve(SigmaABC_inv + ncol(dmx)*Sigma_inv)
    return(list(av = res,
                mu = rmvnorm(1,tmpS%*%apply(Sigma_inv%*%dmx,1,sum),tmpS)))
  },mc.cores = ncores)
  # update random effects muABC
  for (j in 1:nABC){
    av[,hashABC==j] = tmp[[j]]$av
    muABC[,j] = tmp[[j]]$mu
  }
  # update SigmaABC
  SigmaABC = riwish(2*p + nABC, diag(1,p)+muABC%*%t(muABC))
  SigmaABC_inv = solve(SigmaABC)
  if (iter > nb){
    vSigma = cbind(vSigma,c(SigmaABC))
  }
  # Fix mu and Sigma at conditional maximum
  # adaptively update mu and Sigma using MM in the first half of burn-ins
  if (iter <= nb/2){
    # subtract the random effects
    av = mMmu_x(muABC, hashABC, av, FALSE)
    mu = 0.01*apply(av,1,mean)+0.99*mu
    Sigma = 0.01*cov(t(av))+0.99*Sigma
    vmu = cbind(vmu,mu)
    vS = cbind(vS,c(Sigma))
    # add the random effects back
    av = mMmu_x(muABC, hashABC, av, TRUE)
  }
}


# check if the learned parameters make sense
SigmaABC_s = matrix(apply(vSigma,1,mean),3,3)
muABC_s = t(rmvnorm(nABC,rep(0,p),SigmaABC_s))
muxxx = t(rmvnorm(n,mu,Sigma))+muABC_s[,hashABC]
y_s = matrix(rpois(n*p,exp(c(muxxx))),p,n)

# use the wrapper function
res = plgmc(y,hashABC,nABC,100,100,0.99)
