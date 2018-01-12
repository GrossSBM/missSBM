
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
N  <- 5
n  <- 100
Q  <- 2
pi <- diag(.2, Q) +.05
alpha <- rep(1,Q)/Q

# X     <- drawCov(N,n)
# beta  <- drawBeta(N,Q)
# alpha <- drawAlpha(n,Q,beta,X)
###

drawCovSBM_typeII <- function(N,n,Q,pi,X=NULL,beta=NULL,directed=FALSE) {
  if(is.null(X))    X    <- drawCov(N,n)
  if(is.null(beta)) beta <- drawbeta(N,Q)
  alpha <- drawAlpha(n,Q,beta,X)

  Z <- do.call(rbind, lapply(1:n, function(i){ t(rmultinom(1, 1, alpha[i,])) }))
  Znum <- Z %*% c(1:Q)

  Y <- matrix(rbinom(n^2,1,Z %*% pi %*% t(Z)),n)
  if(!directed) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
  diag(X) <- 0

  return(list(Y=Y,cl=as.numeric(Znum),Z=Z, alpha=alpha, beta=beta))
}

### Exemple :
# sbm <- drawCovSBM_typeII(N,n,Q,pi,X,beta)
###

drawTheta <- function(N,m=0,s=1){
  theta  <- rnorm(N,m,s)
  return(theta)
}

sampleCovSBM <- function(Y,X,theta=NULL,beta1=NULL,beta2=NULL,intercept=NULL,directed=FALSE,Node=TRUE){
  n <- nrow(Y)
  if(Node){
    if(is.null(theta)) theta <- drawTheta(nrow(X))
    if(is.null(intercept)) intercept <- rnorm(n)

    sampP  <- 1/(1+exp(-intercept - t(theta) %*% X))
    obsNodes <- which(runif(ncol(X)) < sampP)
    samplingMatrix <- matrix(0,ncol(X),ncol(X)) ; samplingMatrix[obsNodes,] <- 1
    diag(samplingMatrix) <- 1
  } else {
    if(is.null(theta)){
      beta1 <- drawTheta(nrow(X))
      beta2 <- drawTheta(nrow(X))
    }
    if(is.null(intercept)) intercept <- matrix(rnorm(n^2),n,n)

    P <- matrix(0,n,n)
    for(i in 1:n){
      for(j in 1:n){
        P[i,j] <- 1/(1+exp(-intercept[i,j] - t(beta1) %*% X[,i] - t(beta2) %*% X[,j]))
      }
    }
    samplingMatrix <- 1*(P > matrix(runif(n^2),n,n))
    diag(samplingMatrix) <- 1
  }
  if(!directed){
    samplingMatrix <- (t(samplingMatrix) | samplingMatrix)*1
  }
  return(samplingMatrix)
}


### Exemple :
sampMat <- sampleCovSBM(sbm$Y,X,Node=FALSE)
###

drawCovSBM_typeI <- function(N,n,Q,alpha,pi,X=NULL,beta1=NULL,beta2=NULL,directed=FALSE) {
  if(is.null(X))    X    <- drawCov(N,n)
  if(is.null(beta1)){beta1 <- drawBeta(N,Q); beta2 <- drawBeta(N,Q)}

  Z <- t(rmultinom(n, size = 1, prob = alpha))
  Znum <- Z %*% c(1:Q)

  P <- matrix(0,n,n)
  for(i in 1:n){
    for(j in 1:n){
      P[i,j] <- 1/(1+exp(-pi[Znum[i], Znum[j]] - t(beta1) %*% X[,i] - t(beta2) %*% X[,j]))
    }
  }

  Y <- 1*(P > matrix(runif(n^2),n,n))

  if(!directed) Y <- Y * lower.tri(Y) + t(Y * lower.tri(Y))
  diag(X) <- 0

  return(list(Y=Y,cl=as.numeric(Znum),Z=Z, beta1=beta1, beta2=beta2))
}

### Exemple :
# sbm <- drawCovSBM_typeI(N,n,Q,alpha,pi)
###








