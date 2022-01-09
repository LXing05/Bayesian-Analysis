# theta = (theta[1], theta[2])
# target distribution of theta: f(theta) = (2*theta[1]+3*theta[2]+2)/28
# goal: to draw from f(theta)

#################################################################################
#  1) Gibbs Sampler coupled with CDF inversion
#################################################################################
D = 5000
theta = matrix(1, D, 2)

for(i  in 1:(D-1))
{
  # sample theta1 from f(theta1|theta2) using the CDF inversion sampling approach
  u = runif(1,0,1)
  theta[i+1,1]= sqrt(4*u*(6*theta[i,2]+8)+ (3*theta[i,2]+2)^2)/2-(3*theta[i,2]+2)/2
  
  # sample theta2 from f(theta2|theta1) using the CDF inversion sampling approach
  u = runif(1,0,1)
  theta[i+1,2]= sqrt(6*u*(4*theta[i+1,1]+10) + 4*(theta[i+1,1]+1)^2)/3-(2*theta[i+1,1]+2)/3  
}

par(mfrow=c(2,2), mar=c(2,2,1,1))
plot(1:D, theta[,1], type="l")
plot(1:D, theta[,2], type="l")

acf(theta[,1])
acf(theta[,2])

hist(theta[101:D,1], freq=FALSE) # burn-in period 100  
hist(theta[101:D,2], freq=FALSE) # burn-in period 100  


#################################################################################
#  2) Gibbs Sampler coupled with rejection sampling
#################################################################################
D = 5000
theta = matrix(1, D,2)
for(i  in 1:(D-1))
{
  # sample theta1 from f(theta1|theta2) using the rejection sampling approach
  # proposal distribution U(0,2), M=25
  accept=0
  while(accept==0)
  {
    z = runif(1,0,2)
    if( (2*z+3*theta[i,2]+2)>(25*runif(1,0,1)*0.5)){ theta[i+1,1]=z; accept=1}
  }
  
  # sample theta2 from f(theta2|theta2) using the rejection sampling approach
  # proposal dsitribution U(0,2), M=25
  accept=0
  while(accept==0)
  {
    z = runif(1,0,2)
    if( (3*z+2*theta[i+1,1]+2)>(25*runif(1,0,1)*0.5)){ theta[i+1,2]=z; accept=1}
  }
}

par(mfrow=c(2,2), mar=c(2,2,1,1))
plot(1:D, theta[,1], type="l")
plot(1:D, theta[,2], type="l")

acf(theta[,1])
acf(theta[,2])

hist(theta[101:D,1], freq=FALSE) # burn-in period 100  
hist(theta[101:D,2], freq=FALSE) # burn-in period 100 



######################################################################
# example 3:  textbook page 127
# posterior inferenece of the parameters from bivariate normal dist'n
# mu1, mu2, sigma1, sigma2, rho (here sigma1 and sigma2 refer to the variance iso SD)
# prior: 
#  mu1 \propto constant
#  mu2 \propto constant
#  Sigma \propto |Sigma|^{-3/2}
# posterior: 
#   f(mu|Sigma, y) = N(ybar, Sigma/n)
#   f(Sigma|mu,y) = IW(S,v)
#   where v+p+1=n+3, and p=2, so v=n 
######################################################################
# simulate data y
#install.packages("mvtnorm")
library(mvtnorm)
library(MCMCpack)
Sigma = matrix(c(7.5, 0.5*sqrt(7.5*7.0),0.5*sqrt(7.5*7.0), 7.0), 2,2)
n= 200
set.seed(10)
y = round(rmvnorm(n, c(200,150),Sigma),1) 


#######################
n= nrow(y)
ybar = apply(y,2,mean)

M=10000
mu = matrix(100,M,2)  #mean vector
Sigma = array(1,dim=c(M,2,2)); 
Sigma[,1,2]=Sigma[,2,1]=0  # covariance matrix
for(i in 2:M)
{
  #draw mu given Sigma
  mu[i,]= rmvnorm(1, ybar, Sigma[i-1,,]/n)
  
  #draw Sigma given mu
  e = y - matrix(rep(mu[i,],n),byrow=T,ncol=2)
  SS = t(e)%*%(e) # sum of squares error
  Sigma[i,,]<- riwish(n,SS)
}

rho =Sigma[,1,2]/sqrt(Sigma[,1,1]*Sigma[,2,2])

par(mfrow=c(3,2),mar=c(2,2,1,1))
plot(1:M, mu[,1], type="l");
plot(1:M, mu[,2], type="l");
plot(1:M, Sigma[,1,1], type="l");
plot(1:M, Sigma[,2,2], type="l");
plot(1:M, rho, type="l");

par(mfrow=c(3,2),mar=c(2,2,1,1))
acf(mu[1000:M,1])
acf(mu[1000:M,2])
acf(Sigma[1000:M,1,1])
acf(Sigma[1000:M,2,2])
acf(rho[1000:M])

par(mfrow=c(3,2),mar=c(2,2,1,1))
hist(mu[501:M,1]);  
hist(mu[501:M,2]);  
hist(Sigma[501:M,1,1]); 
hist(Sigma[501:M,2,2]);  
hist(rho[501:M]);  
summary(mu[501:M,1]);      quantile(mu[501:M,1], c(0.025,0.975))
summary(mu[501:M,2]);      quantile(mu[501:M,2], c(0.025,0.975))
summary(Sigma[501:M,1,1]); quantile(Sigma[501:M,1,1], c(0.025,0.975))
summary(Sigma[501:M,2,2]); quantile(Sigma[501:M,2,2], c(0.025,0.975))
summary(rho[501:M]);       quantile(rho[501:M], c(0.025,0.975))



