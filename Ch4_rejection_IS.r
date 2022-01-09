#########################################
# Inverse Transformation Sampling
#########################################

# exponential distribution: f(y)=exp(-y/lambda)*lambda
# CDF: 1-exp(-lambda*Y)
# Inverse CDF: y= -ln(1-U)/lambda

lambda=2
u= runif(1000)
y= -log(1-u)/lambda
plot(density(y))
lines(density(rexp(1000,2)), col='red')


# target distribution: f(x)= (2*x+3)/40 for x= [0,5]
# F(x) = (x^2+3x)/40
# F^{-1}(u)  = [-3 + sqrt(9+160p)]/2 since x > -1
u = runif(1000, 0,1)
y = (sqrt(160*u+9)-3)/2
plot(1:1000, y, type="l")
plot(density(y))
summary(y)


#######################################
# rejection Sampling: example 1
# target distribution: f(x)= (2*x+3)/40 
########################################
count=0 # counter for the total number of draws
k=0     # counter for the accepted draws
x = rep(0,5000) # vector to store the accepted draws

while(k<5000)
{
  count = count+1
  z= runif(1,0,5) # g(x)= Unif(0,5)=1/5
  Mgz = 5*(1/5)   # M*g(x), where M=5
  fz = (2*z+3)/40
  
  u= runif(1)
  if(u<fz/Mgz)  
  { 
    k = k +1
    x[k] = z
  }
}
count 
5000/count  # acceptant rate ; should be around 0.2
par(mfrow=c(1,2), mar=c(2,2,1,1))
plot(seq(0,5,0.01), (2*seq(0,5,0.01)+3)/40, type="l")
hist(x)
summary(x)


############################################################################
#  rejection Sampling: example 2
#  target function: defined by function post.dist
############################################################################

# target distribution: f(x)= exp(-a|x-mu|/b)/(2b)
# mean= mu; variane= 2*b^2, support = (-inf, inf)
# envolpe distribution: g(x)= N(mu,b) 
fx = function(x, mu,b) 
{
  target = exp(-abs(x-mu)/b)/(2*b) 
  return(target)
}   
  
# proposal distribution
gx= function(x,mu,b) 
{
  target = exp(-1/2*(x-mu)^2/(2*b^2))/sqrt(2*pi*2*b^2)
  return(target)
}   

# plot the fx and calculate the maximum density
x= seq(-10, 10, 0.01); 
fxx= fx(x,1,2)
gxx= gx(x,1,3)
plot(x, fxx, type="l")
lines(x, gxx, type="l", col='red')

M= max(fxx)/max(gxx)
lines(x, gxx*M, col='blue')


reject = function(x, M, S=5000)
{
  samples=rep(0,S)
  i=0      # counter for the accepted draws
  count=0  # counter for the NO of proposals/trieds
  while(i<S)
  {
    z = rnorm(1,1,sqrt(2*3^2)) # proposal dist N(1,sqrt(2*3^2))
    Mgx = M*gx(z,1,3)
    fx = fx(z,1,2)
    u = runif(1)
    if(u<fx/Mgx) 
    {
      i=i+1
      samples[i]=z 
    }
    count=count+1
  }
  list(acceptrate=S/count,samples=samples)
}
out= reject(x, M) 
out$acceptrate

par(mfrow=c(1,2), mar=c(2,2,1,1))
hist(out$samples)
plot(density(out$samples))
x= seq(-10, 10, 0.01);  fxx= fx(x,1,2); lines(x, fxx, type="l", col='red')


############################################################################
#  importance Sampling
#  goals: 1)obtain E(x); 2) draw from f(x)
############################################################################
# target distribution:
# pdf f(x) propto sqrt((a*x)^a*b^b/(a*x+b)^(a+b))/x, where x>0
a=5
b=4
# it sometimes helps working on the log-scale to increase the numerical stability
# and to reduce the numerical errors for extremely large numbers and near-zero numbers.
log.f= function(x) 1/2*a*log(a*x)+1/2*b*log(b)-1/2*(a+b)*log(a*x+b)-log(x)
x=seq(0.000001,10,0.01)
plot(x, exp(log.f(x)), ylim=c(0, 1.5))

# proposal p(x): log-normal (on the log-scale) 
log.p = function(t,mu, sigma) -log(t)-(log(t)-mu)^2/(2*sigma^2)-log(sigma*sqrt(2*pi))
lines(x, exp(log.p(x,-1,5)), col='red')

## draw 5000 samples from the proposal distrbution log-N(-1,5)
M = 5000
U = rnorm(M, -1,5)
# log importance ratio
log.w = log.f(exp(U)) - log.p(exp(U),-1,5) 
# importance ratio
w = exp(log.w)

# the normalization constant
mean(w)  

# E(x)
mean(w*exp(U))/mean(w)

# SIR
rs = sample(exp(U), 1000, w, replace=FALSE)
hist(rs)
plot(density(rs), col='blue')
# It turns out f(x) is the kernel from the F distribution 
# with degrees of freedom a and b !!
lines(seq(0,50,0.1), df(seq(0,50,0.1),a,b), col='red')












