rm(list=ls())
library(missSBM)

### A SBM model : ###
N <- 300
Q <- 3
M <- 10
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE
covariates <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)

mySBM <- SBM_sampler$new(directed, N, alpha, pi, covariates, covarParam, function(x, y) {-abs(x - y)})
mySBM$rBlocks()
mySBM$rAdjMatrix()
mySBM
