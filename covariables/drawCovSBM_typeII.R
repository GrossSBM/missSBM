
drawCov <- function(N,n,m=0,s=1){
  X <- matrix(rnorm(n*N,m,s), nrow=N, ncol=n)
  return(X)
}
drawBeta <- function(N,Q,m=0,s=1){
  beta  <- matrix(rnorm(N*(Q-1),m,s), nrow=N, ncol=Q-1)
  return(beta)
}
drawAlpha <- function(n,Q,beta,X){
  alpha <- matrix(0,n,Q)
  for(i in 1:n){
    for(q in 1:Q){
      if(q != Q){
        alpha[i,q] <- exp(beta[,q] %*% X[,i])/(1 + sum(exp(t(beta) %*% X[,i])))
      } else {
        alpha[i,q] <- 1/(1 + sum(exp(t(beta) %*% X[,i])))
      }
    }
  }
  return(alpha)
}

### Exemple :
N <- 5
n <- 20
Q <- 2
pi <- diag(.45, Q) +.05

X     <- drawCov(N,n)
beta  <- drawBeta(N,Q)
alpha <- drawAlpha(n,Q,beta,X)
###

drawCovSBM <- function(N,n,Q,pi,X=NULL,beta=NULL,directed=FALSE) {
  if(is.null(X))    X    <- drawCov(N,n)
  if(is.null(beta)) beta <- drawbeta(N,Q)
  alpha <- drawAlpha(n,Q,beta,X)
  
  Z <- do.call(rbind, lapply(1:n, function(i){ t(rmultinom(1, 1, alpha[i,])) }))
  Znum <- Z %*% c(1:Q)
  
  Y <- matrix(rbinom(n^2,1,Z %*% pi %*% t(Z)),n)
  if(!directed) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
  diag(X) <- 0
  
  return(list(Y=Y,cl=as.numeric(Znum),Z=Z))
}

### Exemple :
sbm <- drawCovSBM(N,n,Q,pi,X,beta)
###
