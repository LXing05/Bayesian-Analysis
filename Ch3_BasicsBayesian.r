#################################################################
# example 1:
# inference on mu and sigma for the normal distribution using 
# noninformative prior  f(mu, sigma2) prop to 1/sigma2 (Jeffereys prior)
# marginal posterior distribution in sigma2: inverse-gamma distribution 
# marginal posterior distribution in mu: t distributiuon
#################################################################

############## data
y = rnorm(10,2,2.5)
# statistics that will be used later in the posterior dist'n
s2= var(y)
ybar = mean(y)
n= length(y)

############ posterior inference on sigma2:
# draw sigma2^{-1} from gamma distribution and then take the inverse
sigma2.inv = rgamma(10000, (n-1)/2, rate=(n-1)/2*s2)
sigma2 = 1/sigma2.inv
hist(sigma2)
summary(sigma2)
quantile(sigma2,c(0.025, 0.975))

############# posterior samples on std deviation, sigma
sigma = sqrt(sigma2)
hist(sigma)
summary(sigma)
quantile(sigma,c(0.025, 0.975))

############ posterior samples on mu
# approach 1: draw sigma2 first, then draw mu given sigma2
mu1 = rep(0,10000)
for(i in 1:10000) mu1[i] = rnorm(1,ybar, sqrt(sigma2[i]/n))
hist(mu1)
summary(mu1)
quantile(mu1,c(0.025, 0.975))

# approach 2: draw mu directly from the marginal distriution
mu2 = rt(10000, n-1)*sqrt(s2/n)+ybar


plot(density(mu1))
lines(density(mu2), col='red')

#################################################################
# example 2: Box office example on Dirichlet distribution
# a survey:
# movie A: 125
# movie B: 95
# movie C: 75
# movie D: 70
# movie E: 55
# Multinomial likelihood
# Dirichlet Priors: D(1,1,1,1,1)
# Posterior Priors: D(n1+1, n2+1, n3+1, n4+1, n5+1)
#################################################################
#install.packages("gtools")
library(gtools)

post.p = rdirichlet(2000, c(126,96,76,71,56))
p3 = post.p[,3]
summary(p3); hist(p3); quantile(p3,c(0.025,0.975))

# sum of p1 and p2
p12 = post.p[,1]+post.p[,2]
summary(p12); hist(p12); quantile(p12,c(0.025,0.975))

# difference between p5 and p4
p5_4 = post.p[,5]-post.p[,4]
summary(p5_4); hist(p5_4); quantile(p5_4,c(0.025,0.975))

# how likely p1 > 0.25
p1gt0.25 = (post.p[,1]>0.25)
mean(p1gt0.25)


######################################################################
# example 3:  
# posterior inferenece of the parameters from bivariate normal dist'n
# mu1, mu2, sigma1, sigma2, rho (here sigma1 and sigma2 refer to the variance iso SD)
# prior: 
#  mu1 \propto constant
#  mu2 \propto constant
#  Sigma \propto |Sigma|^{-3/2}
############# #########################################################

# simulate data y
#install.packages("mvtnorm")
library(MCMCpack)
library(mvtnorm)
Sigma = matrix(c(7.5, 0.5*sqrt(7.5*7.0),0.5*sqrt(7.5*7.0), 7.0), 2,2)
n= 200
set.seed(10)

############## data
y = round(rmvnorm(n, c(200,150),Sigma),1) 
n= nrow(y)

# used in the conditional posterior dist'n of mu1 and mu2
ybar = apply(y,2,mean)
SS= cov(y)*(n-1)

# draw posterior samples on my and Sigma
# f(Sigma|y) = IW(n-1, SS^-1)
# f(mu|Sigma, y) = N(ybar, Sigma/n)
M=10000
mu = matrix(100,M,2)  #mean vector
Sigma = array(1,dim=c(M,2,2)); 
for(i in 2:M) {
  Sigma[i,,]<- riwish(n-1,SS)
  mu[i,]= rmvnorm(1, ybar, Sigma[i,,]/n)
}

# posteriro inferences 
par(mfrow=c(3,2),mar=c(2,2,1,1))
hist(mu[,1]);  
hist(mu[,2]);  
hist(Sigma[,1,1]); 
hist(Sigma[,2,2]);  
rho =Sigma[,1,2]/sqrt(Sigma[,1,1]*Sigma[,2,2]) #correlation between the two y's
hist(rho); 

summary(mu[,1]);      quantile(mu[,1], c(0.025,0.975))
summary(mu[,2]);      quantile(mu[,2], c(0.025,0.975))
summary(Sigma[,1,1]); quantile(Sigma[,1,1], c(0.025,0.975))
summary(Sigma[,2,2]); quantile(Sigma[,2,2], c(0.025,0.975))
summary(rho);         quantile(rho, c(0.025,0.975))




