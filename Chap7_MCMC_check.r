
#################### PSRF #################################
# an example from Chap6
# draw posterior inferenece on correlation coeff
# of the standard bivariate normal distribution
# (mu1=mu2=0, sigma1=sigma2=1)
# prior: f(\rho) \propto  (1-rho^2)(-3/2)
##########################################################

library(MASS)

# simulate data y
Sigma = matrix(c(1,0.5,0.5, 1), 2,2)
n= 200
y = round(mvrnorm(n, c(0,0),Sigma),1)


# MH algorithm to draw posterior samples of rho
# log(posterior)
logpost = function(rho,n=200){ 
  -(n+3)/2*log(1-rho^2)-
    0.5*(sum(y[,1]^2)-2*rho*sum(y[,1]*y[,2])+sum(y[,2]^2))/(1-rho^2)}

# posterior samples
MH.rho <- function(init.val= 0, M=1000)
{
  rho = matrix(init.val,M); 
  reject = 0
  for(i in 2:M)
  {
    # proposal distribution: U(rho[i-1]-0.07, rho[i-1]+0.07)
    rho.c = rho[i-1] + runif(1,-0.07,0.07)
    if(abs(rho.c)>1){
      reject = reject+1
      rho[i] = rho[i-1]}
    else{
      u = runif(1,0,1)
      if(logpost(rho.c)-logpost(rho[i-1]) >= log(u) )
        rho[i] = rho.c
      else{
        reject = reject +1
        rho[i] = rho[i-1]}
    }
  }
  print(1-reject/(M-1)) # acceptance rate
  return(rho)
}

# first chain
par(mfrow=c(2,2),mar=c(3,3,1,1))
rho.c1<- MH.rho(-0.2)
plot(1:1000, rho.c1, type="l"); 
acf(rho.c1)

# second chain
rho.c2<- MH.rho(0.5)
plot(1:1000, rho.c2, type="l"); 
acf(rho.c2)

# psrf
library(coda)

mh.list <- mcmc.list(list(mcmc(rho.c1[201:1000]),mcmc(rho.c2[201:1000])))
plot(mh.list); 
summary(mh.list); 
autocorr.plot(mh.list)
gelman.diag(mh.list)
# A potential problem with gelman.diag is that it may mis-diagnose convergence 
# if the shrink factor happens to be close to 1 by chance. 

# By calculating the shrink factor at several points in time, gelman.plot shows 
# if the shrink factorhas really converged, or whether it is still fluctuating.
gelman.plot(mh.list)

