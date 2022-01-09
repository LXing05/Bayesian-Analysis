#################################################################
# use EM to infer membership in the Gaussain Mixture model setting
#################################################################


library(mvtnorm)
set.seed(100)
Sigma1 = matrix(0.1*c(1,0,0,1),2,2)
Sigma2 = matrix(0.1*c(1,1,1,4),2,2)
x1 = rmvnorm(100, c(0,0), Sigma1)
x2 = rmvnorm(50,c(-1,-2), Sigma2)
X= rbind(x1,x2)
plot(X)


#######################################
# the EM algorithm
#######################################

Q = function(n, K, X, Sigma, mu, pi, Tij){
  Q=0
  for(i in 1:n){
    for(j in 1:K){
      Q = Q+ Tij[i,j]*log(pi[j])-
        t(X[i,]-mu[j,])%*%solve(Sigma[,,j])%*%(X[i,]-mu[j,])
    }
  }
  return(Q)
}

GMM.EM = function(params0, X, K = 2, tol=.00001, maxits=100, crit=1)
{       
  # Starting points
  n = nrow(X)
  p= ncol(X)
  mu = params0$mu
  Sigma = params0$Sigma
  pi = params0$pi
  
  # probability that Xi belong to each category (gamma_ij in the slides)
  Tij = matrix(0, n, K)
  it = 0
  converged = FALSE
  
  while ((!converged) & (it < maxits)) { 
    print(it)
    
    ### the E-step (compute the expected values of the latent variables)
    for(i in 1:n){
      for (j in 1:K){
        Tij[i,j] = pi[j] * dmvnorm(X[i,], mu[j,], Sigma[,,j])
      }
    }
    Tij = Tij/rowSums(Tij) # print(Tij)
    
    piOld = pi
    muOld = mu
    SigmaOld = Sigma
    parmlistold = c(piOld, muOld, SigmaOld)
    Qt= Q(n, K, X, Sigma, mu, pi, Tij)
    
    ### the M-step
    Tj = colSums(Tij)
    pi = Tj/n
    for(j in 1:K){
      mu[j,]= t(X) %*% Tij[,j]/Tj[j]      
      S= matrix(rep(0,p*p),p,p)
      for(i in 1:n)  S = S + (X[i,]-mu[j,])%*%t(X[i,]-mu[j,])*Tij[i,j]/Tj[j]
      Sigma[,,j]=S
    }
    parmlistcurrent = c(pi, mu, Sigma)
    
    for(i in 1:n){
      for (j in 1:K){
        Tij[i,j] = pi[j] * dmvnorm(X[i,], mu[j,], Sigma[,,j])
      }
    }
    Tij = Tij/rowSums(Tij)
    Qt1= Q(n, K, X, Sigma, mu, pi, Tij)
    
    it = it + 1
    if(crit==1){
      d = abs(Qt1-Qt)
      converged = (d <= tol)
    }
    else {
      d= max(abs(parmlistold - parmlistcurrent))
      converged = (d <= tol)
    }
  }
  # majority rule, in a 2-cluster case, 50% is the cutoff
  group = which(round(Tij)== max(Tij), arr.ind=T)
  group = group[order(group[,1]), 2]          
  out = list(pi=pi, mu=mu, Sigma=Sigma, Tij=Tij, cluster=group, iter=it, d=d) 
} 
  

params0 = list(mu=matrix(c(0,1,0,-1),2,2), Sigma=array(rep(diag(2),2),c(2,2,2)), pi=c(0.5,0.5))
result= GMM.EM(params0, X)
result

# because the data are simulated, true labels are known, the clustering
# accuracy can be checked, which cannot be done in real life prolems. 
(sum(result$cluster[1:100]==1)+sum(result$cluster[101:150]==2))/150


##################################################
# An alternative: Full Bayesian approach algorithm
# Imputation Step + Posterior Step
###################################################
# posterior distrbution on pi, mu and Sigma
# priors: p(pi) \propto const, 
#         p(mu)\propto constat; 
#         p(Sigma) \propto |Sigma|^{-3/2}

library(MCMCpack)
GMM.IP = function(X, K=2, T=100)
{     
  # Starting points
  K = 2
  n = nrow(X)
  p= ncol(X)
  mu = array(0,c(T,p,K))
  Sigma = array(0,c(p,p,T,K))
  for( j in 1:K ) Sigma[,,1,j] = diag(p)

  pi = matrix(1/p,T,K)
  Z = sample(1:2, n,replace=T)
  Z= do.call(cbind, replicate(T, Z, simplify=FALSE)) 
  
  for(t in 1:(T-1))
  {
    print(t)
    # impute Z from f(z|x,theta)\propo f(z,x,\theta)=f(x|z,\theta)f(z|\theta)
    Tij = matrix(0, n, K)
    for(i in 1:n){
      for (j in 1:K)
        Tij[i,j] = pi[t,j] * dmvnorm(X[i,], mu[t,,j], Sigma[,,t,j])
    }
    Tij = Tij/rowSums(Tij)
    #print(Tij)
    for(i in 1:n) Z[i,t+1]= which(rmultinom(1,1,Tij[i,])==1)
    
    
    # draw mu,Sigma, pi from their posterior dist'n
    m = rep(0,K)
    for(j in 1:K) m[j]= sum(Z[,t+1]==j)
    
    pi[t+1,] = rdirichlet(1,m+1)
    for(j in 1:K){
      e = X[Z[,t+1]==j,]-mu[t,,j]
      SS = t(e)%*%(e)
      Sigma[,,t+1,j]<- riwish(m[j],SS)
      mu[t+1,,j] = rmvnorm(1,apply(X[Z[,t+1]==j,],2,mean),Sigma[,,t+1,j]/m[j])
    }
  }
  out = list(pi=pi, mu=mu, Sigma=Sigma, Z=Z) 
} 

K=2; T= 5000; n= 150
result= GMM.IP(X,K,T)

par(mfrow=c(2,3), mar=c(2,2,1,1))
for(j in 1:K){
  plot(1:T,result$pi[,j], type="l")
  plot(1:T,result$mu[,1,j], type="l")
  plot(1:T,result$mu[,2,j], type="l")
}
for(j in 1:K){
  plot(1:T,result$Sigma[1,1,,j], type="l")
  plot(1:T,result$Sigma[1,2,,j], type="l")
  plot(1:T,result$Sigma[2,2,,j], type="l")
}

for(j in 1:K){
  acf(result$pi[,j])
  acf(result$mu[,1,j])
  acf(result$mu[,2,j])
}
for(j in 1:K){
  acf(result$Sigma[1,1,,j])
  acf(result$Sigma[1,2,,j])
  acf(result$Sigma[2,2,,j])
}
# the acf plots suggest some burning and thinning is needed
thin = seq(51,T,25)
post.pi = result$pi[thin,]
post.mu = result$mu[thin,,]
post.Sigma = result$Sigma[,,thin,]
post.Z = result$Z[,thin]
for(j in 1:K){
  acf(post.pi[,j])
  acf(post.mu[,1,j])
  acf(post.mu[,2,j])
}
for(j in 1:K){
  acf(post.Sigma[1,1,,j])
  acf(post.Sigma[1,2,,j])
  acf(post.Sigma[2,2,,j])
}


plot(1:T,result$Z[1,], type="l")
plot(1:T,result$Z[20,], type="l")
plot(1:T,result$Z[50,], type="l")
plot(1:T,result$Z[80,], type="l")
plot(1:T,result$Z[120,], type="l")
plot(1:T,result$Z[150,], type="l")

cluster= rep(0,n)
for(i in 1:n) 
  cluster[i] = ifelse(mean(result$Z[i,100:T]==1)>0.5,1,2)

(sum(cluster[1:100]==1)+sum(cluster[101:150]==2))/150


#################### VI ######################
# install.packages("LaplacesDemon")
library(LaplacesDemon)

# Demon Data
data(demonsnacks)
# Late one night, after witnessing Laplace's Demon in action, I followed him back to what seemed to be his lair. 
# Minutes later, he left again. I snuck inside and saw something labeled 'Demon Snacks'. 
# Hurriedly, I recorded the 39 items, each with a name and 10 nutritional attributes.
head(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,10]+1)))
J <- ncol(X)
for (j in 2:J) X[,j] <- CenterScale(X[,j])


# Data List Preparation 
# The list of data must include mon.names which contains monitored variable names, 
# parm.names which contains parameter names. 
mon.names <- "mu[1]"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)

MyData <- list(J=J, X=X, mon.names=mon.names,
     parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)

	 
	 
#################  Model Specification  ##########################
Model <- function(parm, Data)
{
     ### Parameters
     beta <- parm[Data$pos.beta]
     sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf) #sigma>0
     parm[Data$pos.sigma] <- sigma
	 
     ### Log-Prior
     beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
     sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
	 
     ### Log-Likelihood
     mu <- tcrossprod(Data$X, t(beta))
     LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
	 
     ### Log-Posterior
     LP <- LL + beta.prior + sigma.prior
	 
     Modelout <- list(LP=LP, Dev=-2*LL, Monitor=mu[1],
          yhat=rnorm(length(mu), mu, sigma), parm=parm)
     return(Modelout)
     }

############################  Initial Values  #############################
Initial.Values <- rep(0,J+1)
# assumes p(theta) ~ N(mean,var) and tries to optimize mean and var
# that maximizes the ELBO
Fit <- VariationalBayes(Model, Initial.Values, Data=MyData, Covar=NULL,
     Iterations=5000, Method="Salimans2", Stop.Tolerance=0.01)

names(Fit)
for(j in 1:(J+1)) hist(Fit$Posterior[,j])
print(Fit)

Pred <- predict(Fit, Model, MyData, CPUs=1)
plot(Pred, Style="Density", Rows=1:9)
plot(Pred, Style="Fitted")



######################################################################

#very slow or just freezes
if(0){
  library(keras)
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()
  K <- keras::backend()
}

library(ruta)
library(purrr)


x <- iris[, 1:4] 
x_train <- x[1:100, ]
x_test <- x[101:150, ]

autoencoder(
  input() + dense(256) + dense(36, "tanh") + dense(256) + output("sigmoid"),
  loss = "mean_squared_error"
) %>%
  make_contractive(weight = 1e-4) %>%
  train(x_train, epochs = 40) %>%
  evaluate_mean_squared_error(x_test)