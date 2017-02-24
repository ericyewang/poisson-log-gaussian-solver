plgmcmc <- function(y,hashABC,nABC,nb,ns,rho){
  cat("loading dependencies...")
  suppressMessages(require(parallel))
  suppressMessages(require(mvtnorm))
  suppressMessages(require(MCMCpack))
  suppressMessages( require(Rcpp) )
  source("~/Dropbox/MyGitbub/poisson-log-gaussian-solver/ars_fcn.R")
  sourceCpp("~/Dropbox/PhD Research/Hierarchical Random Effect Model_Statistical Ecology/HMCEM/monte_carlo_int.cpp")
  cat("finished\n")
  
  # 1. preparation
  cat("preparing for MCMC...")
  n = ncol(y)
  p = nrow(y)
  # 1.1. preparation for parallel computing
  y_listABC = list()
  for (j in 1:nABC){
    y_listABC[[j]] = y[,hashABC==j]
  }
  ncores = detectCores()-1 # number of cores to use
  # 1.2. initial estimate of mu and Sigma
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
  cat("finished\n")
  
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
      mu = (1-rho)*apply(av,1,mean)+rho*mu
      Sigma = (1-rho)*cov(t(av))+rho*Sigma
      vmu = cbind(vmu,mu)
      vS = cbind(vS,c(Sigma))
      # add the random effects back
      av = mMmu_x(muABC, hashABC, av, TRUE)
    }
  }
  return(list(vmu=vmu,vS=vS,vSigma=vSigma,
              mu=mu,Sigma=Sigma))
}