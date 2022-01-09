
#################### MH algorithms ############################
# example 1: 
# Likelihood: f(y_i|theta) \propto theta*y_i+(2-25*theta)/10 for i=1,2,...,n=1350
# where Y={0,1,2,3,4,5} and theta \in (-2/25, 2/25) so that f(y|theta) > 0
# prior f(theta) = U(-2/25, 2/25) \propto const  
# sample from posterior f(\theta|y) \propto likelihood* prior \propto likelihood


# data
y = rep(0,1350) 
y[1:20]=0   
y[21:35]=1   
y[36:100]=2
y[101:425]=3 
y[426:825]=4
y[826:1350]=5


# likelihood
fy = function(y,theta){
  den=1
  for(i in 1:length(y)){
    den= den*(theta*y[i]+(2-25*theta)/10)
  }
  return(den)
}

# M=5000 posterior draws
# proposal/jumping distribution: N(theta[i-1], sigma=0.002)
reject = 0
M = 5000
theta = rep(0.1, M)
for(i in 2:M)
{
  OK=FALSE
  while(!OK) {
    thetac = rnorm(1,theta[i-1],sigma)
    OK=(thetac < 2/25 & thetac > -2/25)
  }
  A = fy(y,thetac)/fy(y,theta[i-1])  # calculate A
  u = runif(1,0,1) 
  if(A>u) theta[i]= thetac
  else{
    theta[i]= theta[i-1]
    reject= reject+1
  }
}
print(1-reject/M)
plot(1:M, theta)



# the above DOES NOT WORK because of computer underflow
# solution: work on the log-scale
# log-likelihood
logfy = function(y,theta){  
  logden= sum(log(theta*y+(2-25*theta)/10))
  return(logden)
}

# parameter m, M=5000 posterior draws
# proposal/jumping distribution: N(theta[i-1], sigma)
eg1 <- function(sigma, M=5000, thin=10, burnin=500){
  reject = 0
  theta = rep(0, M) # m[1] = 0

  for(i in 2:M)
  {
    OK=FALSE
    while(!OK) {
      thetac = rnorm(1,theta[i-1],sigma)
      OK=(thetac < 2/25 & thetac > -2/25)
    }
    logA = logfy(y,thetac)-logfy(y,theta[i-1])
    logu = log(runif(1,0,1)) 
    if(logA>logu)
      theta[i]= thetac
    else{
      theta[i]= theta[i-1]
      reject= reject+1}
    
  }
  par(mfrow=c(2,2),mar=c(2,2,1,1))
  plot(1:M, theta, type="l");
  theta = theta[burnin:M]
  
  plot(acf(theta, plot=FALSE))
  pick = seq(1,M-burnin,thin)
  
  theta = theta[pick]
  plot(acf(theta, plot=FALSE))
  hist(theta);    
  print(summary(theta))
  return(list(theta=theta, A = 1-reject/(M-1)))
}
# acceptance rate:
draw2<- eg1(0.002); draw2$A;           
draw1<- eg1(0.0005, thin=30); draw1$A; 
draw3<- eg1(0.1, thin=35); draw3$A;    



######################################################################
# example 2: 
# draw posterior inferenece on correlation coeff
# of the standard bivariate normal distribution
# (mu1=mu2=0, sigma1=sigma2=1 known)
# prior: f(\rho) \propto  (1-rho^2)^(-3/2)
######################################################################

library(MASS)

# simulate data y
Sigma = matrix(c(1,0.5,0.5, 1), 2,2)
n= 200
y = round(mvrnorm(n, c(0,0),Sigma),1)

# log(posterior)
logpost = function(rho,n=200){ 
  -(n+3)/2*log(1-rho^2)-
    0.5*(sum(y[,1]^2)-2*rho*sum(y[,1]*y[,2])+sum(y[,2]^2))/(1-rho^2)}

# posterior samples
eg2 <- function(width, M=5000){
  rho = matrix(0,M); 
  reject = 0
  # this is the tuning parameter for jumping width
  for(i in 2:M){
    # proposal distribution: U(rho[i-1]-0.07, rho[i-1]+0.07)
    OK=FALSE
    while(!OK) {
      rho.c = rho[i-1] + runif(1,-width,width)
      OK= (abs(rho.c)<=1)
    }
    u = runif(1,0,1)
    if(logpost(rho.c)-logpost(rho[i-1]) >= log(u) )
      rho[i] = rho.c
    else{
      reject = reject +1
      rho[i] = rho[i-1]}
    
  }
  return(list(rho=rho, A = 1-reject/(M-1)))
}

draw1<- eg2(0.07); draw1$A
par(mfrow=c(2,2),mar=c(2,2,1,1))
plot(1:M, draw1$rho, type="l");
acf(draw1$rho)

pick = seq(1,4500,20)
rho = draw1$rho[501:5000][pick]
acf(rho); hist(rho);  
summary(rho); quantile(rho, c(0.025,0.975))



#  Gibbs Sampler
#################################################################################
# example 1: GS coupld with inverse-transformation sampling
# theta = (theta[1], theta[2]) \in [0,2]x[0,2]
# target distribution of theta: f(theta) = (2*theta[1]+3*theta[2]+2)/28
# goal: to draw from f(theta) 
#################################################################################
D = 5000
theta = matrix(1, D, 2)
t0= proc.time()
for(i  in 1:(D-1))
{
  # sample theta1 from f(theta1|theta2) using the CDF inversion sampling approach
  u = runif(1,0,1)
  theta[i+1,1]= sqrt(4*u*(6*theta[i,2]+8)+ (3*theta[i,2]+2)^2)/2-(3*theta[i,2]+2)/2
  
  # sample theta2 from f(theta2|theta1) using the CDF inversion sampling approach
  u = runif(1,0,1)
  theta[i+1,2]= sqrt(6*u*(4*theta[i+1,1]+10) + 4*(theta[i+1,1]+1)^2)/3-(2*theta[i+1,1]+2)/3  
}
proc.time()-t0
par(mfrow=c(2,2), mar=c(2,2,1,1))
plot(1:D, theta[,1], type="l")
plot(1:D, theta[,2], type="l")

acf(theta[,1])
acf(theta[,2])

hist(theta[101:D,1], freq=FALSE) # burn-in period 100  
hist(theta[101:D,2], freq=FALSE) # burn-in period 100  

plot(theta[101:D,])
cor(theta[101:D,])


#################################################################################
#  example 2: Gibbs Sampler coupled with rejection sampling
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


################################################################################
# example 3: 
# Gibbs Sampler drawing mu and Sigma alternatively in Gaussian likelihood.
# posterior inferenece of the parameters from bivariate normal dist'n
# mu1, mu2, sigma1, sigma2, rho (here sigma1 and sigma2 refer to the variance iso SD)
# prior: 
#  mu1 \propto constant
#  mu2 \propto constant
#  Sigma \propto |Sigma|^{-3/2} (IW(v=0))
# posterior: 
#   f(mu|Sigma, y) = N(ybar, Sigma/n)
#   f(Sigma|mu,y) = IW(S,v)
#   where v+p+1=n+3, and p=2, so v=n 
################################################################################

# simulate data y
#install.packages("mvtnorm")
library(MASS)
library(MCMCpack)
library(mvtnorm)
Sigma = matrix(c(7.5, 0.5*sqrt(7.5*7.0),0.5*sqrt(7.5*7.0), 7.0), 2,2)
n= 200
set.seed(10)
y = round(mvrnorm(n, c(200,150),Sigma),1) 


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
hist(Sigma[501:M,]);  
hist(rho[501:M]);  
summary(mu[501:M,1]);      quantile(mu[501:M,1], c(0.025,0.975))
summary(mu[501:M,2]);      quantile(mu[501:M,2], c(0.025,0.975))
summary(Sigma[501:M,1,1]); quantile(Sigma[501:M,1,1], c(0.025,0.975))
summary(Sigma[501:M,2,2]); quantile(Sigma[501:M,2,2], c(0.025,0.975))
summary(rho[501:M]);       quantile(rho[501:M], c(0.025,0.975))


###############
# MH + GS
# example 1:  
# posterior inferenece of the parameters from bivariate normal dist'n
# mu1, mu2, sigma1, sigma2, rho (here sigma1 and sigma2 are the variances iso SD)
# prior: 
#  mu1 \propto constant
#  mu2 \propto constant
#  Sigma \propto |Sigma|^{-3/2} (nu-->0 and |Psi|--> 0)
############# #########################################################

# simulate data y
#install.packages("mvtnorm")
library(MCMCpack)
library(MASS)
Sigma = matrix(c(7.5, 0.5*sqrt(7.5*7.0),0.5*sqrt(7.5*7.0), 7.0), 2,2)
n= 200
set.seed(10)
y = round(mvrnorm(n, c(200,150),Sigma),1) 

#################################################
n= nrow(y)
# used in the conditional posterior dist'n of mu1 and mu2
ybar1 = mean(y[,1]); ybar2 = mean(y[,2])


# joint posterior dist'n of sigma1, sigma2, rho, mu1, mu2 on the log scale
logpost = function(ss1, ss2, ss12, rho, sigma1, sigma2, n=200){
  -(3+n)/2*log((1-rho^2)*sigma1*sigma2)-
    0.5/(1-rho^2)*(ss1/sigma1-2*rho*ss12/sqrt(sigma1*sigma2)+ss2/sigma2) 
}
# ss1=sum((y[,1]-mu1[i])^2); 
# ss2=sum((y[,2]-mu2[i])^2); 
# ss12=sum((y[,1]-mu1[i])*(y[,2]-mu2[i]))


logIW = function(X, V, df)
{
  p=nrow(X)
  logg = df/2*log(det(V))-(df+p+1)/2*log(det(X))-sum(diag(V%*%solve(X)))/2
  return(logg)
}

# the outer framework is GS
# When drawing Sigma fromt its full conditional distribution, MH is used
# where the proposal dist'n is IW, and df of the IW is used to tune the acceptance rate.
# more explanation is provided below.
eg1 <- function(df=10, M=10000, burnin=1000, thin=15){
  mu1 = rep(ybar1,M); # starting value at ybar1
  mu2 = rep(ybar2,M); # starting value at ybar2
  Sigma = array(cov(y), c(2,2,M)) # starting value is the sample variance/covariance matrix
  sigma1 = Sigma[1,1,];  
  sigma2 = Sigma[2,2,];  
  rho = Sigma[1,2,]/sqrt(sigma1*sigma2)
  
  reject = 0
  for(i in 2:M){
    
    #sample mu1 from its full conditional dist'n (direct sampling)
    mu1[i]=rnorm(1, ybar1+sqrt(sigma1[i-1]/sigma2[i-1])*rho[i-1]*(mu2[i-1]-ybar2),
                 sqrt(sigma1[i-1]*(1-rho[i-1]^2)/n))
    
    #sample mu2 from its full conditional dist'nl (direct sample)
    mu2[i]=rnorm(1, ybar2+sqrt(sigma2[i-1]/sigma1[i-1])*rho[i-1]*(mu1[i]-ybar1),
                 sqrt(sigma2[i-1]*(1-rho[i-1]^2)/n))
    
    #update sum of squares/cross-product
    ss1=sum((y[,1]-mu1[i])^2); 
    ss2=sum((y[,2]-mu2[i])^2); 
    ss12=sum((y[,1]-mu1[i])*(y[,2]-mu2[i]))
    
    # sample Sigma using (MH)
    # Sigma|mu1, mu2, y \propto f(Sigma, mu1, mu2, rho|y)
    # Let x~ IW(df,S), E(x)=S/(df-1-p)
    # here p =2,  S = E(x)*(df-3)
    # The variance of x is available for each individaul entry, 
    # the large df is, the smaller the variances are
    Sigma.c = riwish(df, (df-3)*Sigma[,,i-1]) # jumping distribution
    #print(Sigma)
    
    u= runif(1,0,1)
    logA1 = logpost(ss1,ss2,ss12,Sigma.c[1,2]/sqrt(Sigma.c[1,1]*Sigma.c[2,2]),Sigma.c[1,1],Sigma.c[2,2])-
      logpost(ss1,ss2,ss12,Sigma[1,2,i-1]/sqrt(Sigma[1,1,i-1]*Sigma[2,2,i-1]),Sigma[1,1,i-1], Sigma[2,2,i-1])
    logA2 =logIW(Sigma[,,i-1], Sigma.c, df)-logIW(Sigma.c, Sigma[,,i-1], df) 
    #print(c(logA1, logA2))
    
    if(logA1+logA2 >log(u))
      Sigma[,,i] = Sigma.c
    else{
      reject = reject+1; 
      Sigma[,,i] = Sigma[,,i-1]}
    sigma1[i] = Sigma[1,1,i]
    sigma2[i] = Sigma[2,2,i]
    rho[i] = Sigma[1,2,i]/sqrt(Sigma[1,1,i]*Sigma[2,2,i])
  }
  
  A= 1-reject/(M-1)
  print(A)
  par(mfrow=c(3,2),mar=c(2,2,1,1))
  plot(1:M, mu1, type="l");
  plot(1:M, mu2, type="l");
  plot(1:M, sigma1, type="l");
  plot(1:M, sigma2, type="l");
  plot(1:M, rho, type="l");
  
  par(mfrow=c(3,2),mar=c(2,2,1,1))
  acf(mu1[(burnin+1):M])
  acf(mu2[(burnin+1):M])
  acf(sigma1[(burnin+1):M])
  acf(sigma2[(burnin+1):M])
  acf(rho[(burnin+1):M])
  
  
  # use the thinning period across all parameters
  pick = seq(1,M-burnin,thin)
  new.mu1 = mu1[(burnin+1):M][pick];
  new.mu2 = mu2[(burnin+1):M][pick]; 
  new.s1 = sigma1[(burnin+1):M][pick];
  new.s2 = sigma2[(burnin+1):M][pick];
  new.rho = rho[(burnin+1):M][pick];
  
  par(mfrow=c(3,2),mar=c(2,2,1,1))
  acf(new.mu1)
  acf(new.mu2)
  acf(new.s1)
  acf(new.s2)
  acf(new.rho)
  
  par(mfrow=c(3,2),mar=c(2,2,1,1))
  hist(new.mu1);  
  hist(new.mu2);  
  hist(new.s1); 
  hist(new.s2);  
  hist(new.rho);  
  
  print(summary(new.mu1)); print(quantile(new.mu1, c(0.025,0.975)))
  print(summary(new.mu2)); print(quantile(new.mu2, c(0.025,0.975)))
  print(summary(new.s1));  print(quantile(new.s1, c(0.025,0.975)))
  print(summary(new.s2));  print(quantile(new.s2, c(0.025,0.975)))
  print(summary(new.rho)); print(quantile(new.rho, c(0.025,0.975)))
}
# eg1(10)
# eg1(100)
# eg1(50)



###############################################################
# Example 2: alternative approach to eg3
# draw mu1 from f(mu1|mu2, Sigma), which is univariate normal
# draw mu2 from f(mu2|mu1, Sigma), which is univariate normal
# draw sigma1^2 f(sigma1|mu1, mu2, sigma2, rho) via MH
# draw sigma2^2 f(sigma1|mu1, mu2, sigma1, rho) via MH
# draw rho from f(rho|mu1, mu2, sigma1, sigma2) via MH
############# #########################################################

n= nrow(y)
# used in the conditional posterior dist'n of mu1 and mu2
ybar1 = mean(y[,1]) ; ybar2 = mean(y[,2])


# w1 and w2 are the jumping widths of two jumping dist'n respectively.
eg2 <- function(w1, w2, M=10000){
  mu1 = matrix(0,M); # starting value 0
  mu2 = matrix(0,M); # starting value 0
  sigma1 = matrix(1,M);  # starting value 1
  sigma2 = matrix(1,M);  # starting value 1
  rho = matrix(0,M) # starting value 0
  reject1 = 0
  reject2 = 0
  reject3 = 0
  for(i in 2:M){
    
    #sample mu1 from its full conditional dist'n (direct sampling)
    mu1[i]=rnorm(1, ybar1+sqrt(sigma1[i-1]/sigma2[i-1])*rho[i-1]*(mu2[i-1]-ybar2),
                 sqrt(sigma1[i-1]*(1-rho[i-1]^2)/n))
    
    #sample mu2 from its full conditional dist'nl (direct sampling)
    mu2[i]=rnorm(1, ybar2+sqrt(sigma2[i-1]/sigma1[i-1])*rho[i-1]*(mu1[i]-ybar1),
                 sqrt(sigma2[i-1]*(1-rho[i-1]^2)/n))
    
    #update sums of squares/cross-product
    ss1=sum((y[,1]-mu1[i])^2); 
    ss2=sum((y[,2]-mu2[i])^2); 
    ss12=sum((y[,1]-mu1[i])*(y[,2]-mu2[i]))
    
    # sample sigma1 (MH)
    # f(sigma1|sigma2,mu1, mu2,rho) \propto f(Sigma, mu1, mu2, rho|y)
    # jumping: U(sigma1[i-1]-2.5,sigma1[i-1]+2.5)
    OK = FALSE
    while(!OK) {
      sigma1.c=sigma1[i-1]+runif(1,-w1,w1)
      OK=(sigma1.c>0)
    }
    u= runif(1,0,1)
    if(logpost(ss1,ss2,ss12, rho[i-1], sigma1.c,  sigma2[i-1])-
       logpost(ss1,ss2,ss12, rho[i-1], sigma1[i-1],sigma2[i-1])
       >log(u))
      sigma1[i] = sigma1.c
    else{
      reject1 = reject1+1; 
      sigma1[i] = sigma1[i-1]}
    
    # sample sigma1 (MH)
    # f(sigma1|sigma2,mu1, mu2,rho) \propto f(Sigma, mu1, mu2, rho|y)
    OK = FALSE
    while(!OK) {
      sigma2.c=sigma2[i-1]+runif(1,-w1,w1)
      OK=(sigma2.c>0)
    }
    u= runif(1,0,1)
    if(logpost(ss1,ss2,ss12,rho[i-1],sigma1[i],sigma2.c)-
       logpost(ss1,ss2,ss12,rho[i-1],sigma1[i],sigma2[i-1])
       >log(u))
      sigma2[i] = sigma2.c
    else{
      reject2 = reject2+1; 
      sigma2[i] = sigma2[i-1]}
  
    # sample rho (MH)
    # f(rho|sigma1,sigma2,mu1, mu2) \propto f(Sigma, mu1, mu2, rho|y)
    OK = FALSE
    while(!OK) {
      rho.c=rho[i-1]+runif(1,-w2,w2) # jumping: U(rho[i-1]-0.15,rho[i-1]+0.15)
      OK=(abs(rho.c)<1)
    }
    u= runif(1,0,1)
    if(logpost(ss1,ss2,ss12,rho.c,sigma1[i],sigma2[i])-
       logpost(ss1,ss2,ss12,rho[i-1],sigma1[i],sigma2[i])
       >log(u))
      rho[i] = rho.c
    else{
      reject3 = reject3+1; 
      rho[i] = rho[i-1]}
  }
  A= c(1-reject1/(M-1), 1-reject2/(M-1), 1-reject3/(M-1))
  print(A)
  
  par(mfrow=c(3,2),mar=c(2,2,1,1))
  plot(1:M, mu1, type="l");
  plot(1:M, mu2, type="l");
  plot(1:M, sigma1, type="l");
  plot(1:M, sigma2, type="l");
  plot(1:M, rho, type="l");
  
  par(mfrow=c(3,2),mar=c(2,2,1,1))
  acf(mu1[1001:M])
  acf(mu2[1001:M])
  acf(sigma1[1001:M])
  acf(sigma2[1001:M])
  acf(rho[1001:M])
  
  
  # use the thinning period across all parameters
  pick = seq(1,M-1000,10)
  new.mu1 = mu1[1001:M][pick];
  new.mu2 = mu2[1001:M][pick]; 
  new.s1 = sigma1[1001:M][pick];
  new.s2 = sigma2[1001:M][pick];
  new.rho = rho[1001:M][pick];
  
  par(mfrow=c(3,2),mar=c(2,2,1,1))
  hist(new.mu1);  
  hist(new.mu2);  
  hist(new.s1); 
  hist(new.s2);  
  hist(new.rho);  
  
  print(summary(new.mu1)); print(quantile(new.mu1, c(0.025,0.975)))
  print(summary(new.mu2)); print(quantile(new.mu2, c(0.025,0.975)))
  print(summary(new.s1));  print(quantile(new.s1, c(0.025,0.975)))
  print(summary(new.s2));  print(quantile(new.s2, c(0.025,0.975)))
  print(summary(new.rho)); print(quantile(new.rho, c(0.025,0.975)))
}
eg2(2.5, 0.1)




