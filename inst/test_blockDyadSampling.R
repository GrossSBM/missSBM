rm(list = ls())
library(missSBM)
library(mclust)

n <- 200
Q <- 3
alpha <- rep(1,Q)/Q
pi <- diag(.45,Q) + .05
directed <- FALSE

sbm <- simulateSBM(n, alpha, pi, directed)
adjacencyMatrix <- sbm$adjacencyMatrix
samplingParameters <- diag(.3,Q) + .3
sampling <- "block_dyad"
sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters, clusters = sbm$memberships)
R <- sampledAdjMatrix$samplingMatrix*1

# Sampling rate
(sum(sampledAdjMatrix$samplingMatrix)-n)/(n*(n-1))

# Inference
vBlocks <- 1:5
infer <- inferSBM(adjacencyMatrix, vBlocks, sampling)

# ICL
ICL <- sapply(infer$models, function(x) x$vICL)
plot(ICL)

# Error
sum((infer$models[[3]]$fittedSampling$parameters - samplingParameters)^2)/sum((infer$models[[3]]$fittedSampling$parameters)^2)
adjustedRandIndex(sbm$memberships, infer$models[[3]]$fittedSBM$memberships)

logB <- function(x,a){
  zero <- .Machine$double.eps
  if(a<zero) a <- zero
  if(a>zero) a <- 1-zero
  return(x*log(a) + (1-x)$log(1-a))
}

iclCalculation <- function(X, R, infer, vBlocks){
  K <- length(vBlocks)
  n <- nrow(X)
  icl <- vector("numeric", length = K)
  for (i in 1:K) {
    loglik <- 0
    pen    <- K*(K+1)*log(n*(n-1)/2) + (K-1)*log(n)
    Tau    <- infer$models[[K]]$fittedSBM$blocks
    for (i in 1:n) {
      for (j in 1:n) {
        if(!is.na(X[i,j])){
          for (q in 1:K) {
            for (l in 1:K) {
              loglik <- loglik + Tau[i,q]*Tau[j,l]( logB(X[i,j],infer$models[[K]]$fittedSBM$connectParam[q,l]) +
                                                      logB(R[i,j],infer$models[[K]]$infer$models[[3]]$fittedSampling$parameters[q,l]) )
            }
          }
        }
      }
    }
  }
}
