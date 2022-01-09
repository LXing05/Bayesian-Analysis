#########################################
# Inverse Transformation Sampling
#########################################

# example 1
# exponential distribution: f(y)=exp(-y*lambda)*lambda
# CDF: 1-exp(-lambda*Y)
# Inverse CDF: y= -ln(1-U)/lambda
lambda=2
p= runif(100000)
y= -log(1-p)/lambda
plot(density(y), lwd=2)
lines(density(rexp(1000000,lambda)), col='red') # true density


# example 2
# target distribution: f(x)= (2*x+3)/40 for x= [0,5]
# F(x) = (x^2+3x)/40
# F^{-1}(p)  = [-3 + sqrt(9+160p)]/2 since x > -1
p = runif(100000, 0,1)
y = (sqrt(160*p+9)-3)/2
hist(y, nclass=50, freq = FALSE)  # empirical density based on the drawn samples
segments(0, 3/40, 5, 13/40, col='red') # true density
plot(density(y), lwd=2)
segments(0, 3/40, 5, 13/40, col='red') # true density



#######################################
# rejection Sampling: example 1
# target distribution: f(x)= (2*x+3)/40 
########################################
A=0   #  the total number of samples
t=1   #  index of the accepted samples for a target of T=5000
x = rep(0,5000) # vector to store the accepted draws

while(t<=5000)
{
  A = A+1
  z= runif(1,0,5) # g(x)= Unif(0,5)=1/5
  Mgx = 5*(1/5)   # M*g(x), where M=5
  fx = (2*z+3)/40
  # accept z with probabilty  = fx/Mgx
  u= runif(1)
  if(u<fx/Mgx){ 
    x[t] = z
    t = t +1}
}
A 
5000/A  # acceptant rate ; should be around 0.2 =1/M
par(mfrow=c(1,2), mar=c(2,2,1,1))
plot(seq(0,5,0.01), (2*seq(0,5,0.01)+3)/40, type="l")
hist(x, freq = FALSE)  # empirical density based on the drawn samples
summary(x)


############################################################################
#  rejection Sampling: example 2
#  target function: defined by function post.dist
############################################################################

# target distribution: f(x)= exp(-|x-mu|/b)/(2b)
# mean= mu; variane= 2*b^2, support = (-inf, inf)
# envolpe distribution: g(x)= N(mu,b) 
fx = function(x, mu,b)  target = exp(-abs(x-mu)/b)/(2*b) 

# two options for the proposal distribution g(x): Gaussian ("N") and Uniform
reject = function(M, S=5000, prop="N")
{
  samples=rep(0,S)
  t=1      # Aer for the accepted draws
  A=0  # Aer for the NO of proposals/trieds
  while(t<=S)
  {
    A=A+1
    if(prop=="N"){
      z = rnorm(1,1,sqrt(2*3^2)) # proposal dist N(1,sqrt(2*3^2))
      Mgx = M*gx(z,1,3)}
    else{
      z= runif(1, -20,20); Mgx=M*(1/40)
    }
    fx = fx(z,1,2)
    u = runif(1)
    if(u<fx/Mgx) {
       samples[t]=z 
       t=t+1}
  }
  list(acceptrate=S/A,samples=samples)
}

############## if the proposal dist' is Gaussian ##############

# calculate M
gx= function(x,mu,sigma)  
  target = exp(-1/2*(x-mu)^2/(2*sigma^2))/sqrt(2*pi*2*sigma^2)
# plot the fx and calculate the maximum density
x= seq(-20, 20, 0.01); 
fxx= fx(x,1,2)
gxx= gx(x,1,3)
plot(x, fxx, type="l")
lines(x, gxx, type="l", col='red')
M= max(fxx)/max(gxx)
lines(x, gxx*M, col='blue')

# perform rejection sampling
t0<-proc.time()
out= reject(M) 
proc.time()-t0
out$acceptrate; 1/M


plot(density(out$samples), lwd=2)
x= seq(-20, 20, 0.01);  fxx= fx(x,1,2); lines(x, fxx, type="l", col='red')
summary(out$samples)
var(out$samples)


############## if the proposal dist' is Unif(-20, 20) ##############
x= seq(-20, 20, 0.01); 
fxx= fx(x,1,2)
gxx= rep(1/40,length(x))
plot(x, fxx, type="l")
lines(x, gxx, type="l", col='red')

M= max(fxx)/max(gxx)
lines(x, gxx*M, col='blue')

t0<-proc.time()
out= reject(M, prop="U") 
proc.time()-t0
out$acceptrate; 1/M

plot(density(out$samples), lwd=2)
x= seq(-20, 20, 0.01);  fxx= fx(x,1,2); lines(x, fxx, type="l", col='red')
summary(out$samples)
var(out$samples)


############################################################################
#  adaptive rejection Sampling: 
############################################################################
library(ars)

# log(f(x)), where f(x) is prop to the target dist'n 
logf<-function(x,mu=0,sigma=1){-1/(2*sigma^2)*(x-mu)^2}
# derivative of log(f(x)) w.r.t x
logfderiv<-function(x,mu=0,sigma=1){-1/sigma^2*(x-mu)}

t0<-proc.time()
n=5000
# x: starings points for which log(f) is claculated, make sure there are points 
#    to both the left and right of the mode of f(x)
# m: # of starting pints (length(x))
# mu and sigma are arguments to function logf
mysample<- ars(n,logf,logfderiv, x=c(-8,0,10), m=3, mu=2,sigma=3)
proc.time()-t0
summary(mysample)
var(mysample)

plot(density(mysample), lwd=2)
x= seq(-20, 20, 0.01); lines(x, dnorm(x,2,3), type="l", col='red')




############################################################################
# importance Sampling + Resampling
# target distribution: f(x) propto sqrt((a*x)^a*b^b/(a*x+b)^(a+b))/x, where x>0
# goals: 1) obtain E(x); 2) draw from f(x)
############################################################################
a=5
b=4
# it sometimes helps working on the log-scale to increase the numerical stability
# and reduce the numerical errors for extremely large numbers and near-zero numbers.
log.f= function(x) 1/2*a*log(a*x)+1/2*b*log(b)-1/2*(a+b)*log(a*x+b)-log(x)
x=seq(0.000001,10,0.01)
plot(x, exp(log.f(x)), ylim=c(0, 0.5), type="l")

# proposal p(x): log-normal 
log.p = function(t,mu, sigma) -log(t)-(log(t)-mu)^2/(2*sigma^2)-log(sigma*sqrt(2*pi))
lines(x, exp(log.p(x,-2,5)), col='red')


## draw 5000 samples from the proposal distrbution logN(-1,5)
M = 5000
U = exp(rnorm(M, -2,5))

# log importance ratio
log.w = log.f(U) - log.p(U,-2,5) 


# importance ratio
w = exp(log.w)

# the normalization constant
mean(w)  

# E(x): use the same set  of U to evaluate the numerator and demoninator
mean(w*U)/mean(w)

# ESS for estimating the mean
# normalize weights
w.n<- w/sum(w)
(ESS2<-1/sum(w.n^2) )


# IS/R with replacement
rs = sample(U, 1000, w, replace = TRUE)
plot(density(rs),lwd=2)

# f(x) is actually the kernel from the F distribution 
# with degrees of freedom a=5 and b=4!! E(x)=df2/(df2-2) = 4/(4-2)=2
lines(seq(0,50,0.1), df(seq(0,50,0.1),a,b), col='red')




