
### Author: Yuanyao Tan (tyuanyao@gmail.com) and Wei Liang (lwlx2h@gmail.com)
### Date: 2020-10-03
### Version 2.0
### Reference: Tan, Y., Wen, X., Lang, W. and Yan, Y. (2020), Empirical Likelihood Weighted Estimation 
### of Avearge Treatment Effects in Randomized Clinical Trials, arXiv preprint arXiv:2008.12989.

# Maximizer over lambda-----------
l <- function(lambda0,U){
  return(sum(log(1+lambda0%*%U)))}

#Algorithm ---------

Newt.p <- function(A){ # A is the matrix of the weights in the constraint, i.e., A^Tp=0.
  
  U <- t(A)
  d <- dim(U)[1]
  N <- ncol(U)
  w <- rep(0,N)
  tag <- 0
  K <- 1000
  
  # # convex hull problem-add two points
  # s <- 1  
  # Ubar <- rowMeans(U)
  # u <- Ubar/norm(Ubar,type="2")
  # cu <- (u%*%var(t(U))%*%u)^(-1/2)
  # u1 <- -c(s*cu)*u
  # u2 <- 2*Ubar+c(s*cu)*u
  # U <- cbind(U,u1,u2)
  
  
  lambda0 <- rep(0,d)
  l0 = l(lambda0,U)
  eps=1e-8  
  gamma0 <- 1  
  
  for(j in 1:K){
    g <-  t((1/(1+lambda0%*%U))%*%t(U))
    v <- as.vector(1/(1+lambda0%*%(U)))^2
    V <- diag(v,N,N)
    dg <- U%*%V%*%t(U) 
    
    e <- -solve(dg)%*%g
    
    if(sqrt(t(e)%*%e)<eps){
      tag <- 1
      lambda <- lambda0
      break;
    }
    
    l1 <- l(as.vector(lambda0),U)
    neg <- sum((1+t(lambda0-gamma0*e)%*%U)<=0)
    
      while(neg>0||l1<l0){
      gamma0 <- gamma0/2
      if(gamma0==0)break
      if(neg==0){l1 <- l(as.vector(lambda0-gamma0*e),U)}
      else{neg <- sum((1+t(lambda0-gamma0*e)%*%U)<=0)}
    }
    
    l0 <- l1
    
    lambda0 <- as.vector(lambda0 - gamma0*e)
    gamma0 <- (j+1)^(-1/2)
  }
  
  if(tag==0){
    cat("fail to converge")
    lambda=lambda0
  }
  
  # weight (empirical probability)
  
  w <- as.vector(1/(N*(1+lambda%*%U)))
  return(list(w=w,lambda=lambda,tag=tag,U=t(U)))
}

