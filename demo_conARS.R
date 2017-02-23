# Demo

# Simulation----
y = c(0,0,0)
mu = c(-3,-2,-4)
p = 3 # dimension
Fc = 0.8*matrix(rnorm(p*p),p,p) # loadings for the covariance matrix
Sigma = Fc%*%t(Fc)+diag(0.1,p) # create a covariance matrix
res1 = imh(y,mu,Sigma,1000)
res2 = pm(y,mu,Sigma,1000)
res3 = conARS(y,mu,Sigma,1000)
plot(1:1000,res1$samples[3,],'l')
plot(1:1000,res2$samples[3,],'l')
plot(1:1000,res3[3,],'l')