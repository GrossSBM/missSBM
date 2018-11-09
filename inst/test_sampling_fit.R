rm(list=ls())
library(missSBM)

### A SBM model : ###
n <- 400
Q <- 5
alpha <- rep(1,Q)/Q                       # mixture parameter
pi <- diag(.45,Q) + .05                   # connectivity matrix
directed <- FALSE                         # if the network is directed or not

### Draw a SBM undirected model
mySBM <- simulateSBM(n, alpha, pi, directed)

adjacencyMatrix <- mySBM$adjacencyMatrix

## Draw
sampledNet_dyad <- samplingSBM(adjacencyMatrix, "dyad", 0.1)
fittedSampling_dyad <- missSBM:::dyadSampling_fit$new(sampledNet_dyad)
fittedSampling_dyad$vExpec
fittedSampling_dyad$parameters
fittedSampling_dyad$df
fittedSampling_dyad$penalty
-2 * fittedSampling_dyad$logLik + fittedSampling_dyad$penalty

sampledNet_node <- samplingSBM(adjacencyMatrix, "node", 0.1)
fittedSampling_node <- nodeSampling_fit$new(sampledNet_node)
fittedSampling_node$vExpec
fittedSampling_node$parameters
fittedSampling_node$df
fittedSampling_node$penalty
-2 * fittedSampling_node$logLik + fittedSampling_node$penalty

sampledNet_double_standard <- samplingSBM(adjacencyMatrix, "double_standard", c(0.1, 0.5))
fittedSampling_double_standard <- doubleStandardSampling_fit$new(sampledNet_double_standard)
fittedSampling_double_standard$vExpec
fittedSampling_double_standard$parameters
fittedSampling_double_standard$df
fittedSampling_double_standard$penalty
-2 * fittedSampling_double_standard$logLik + fittedSampling_double_standard$penalty

sampledNet_block <- samplingSBM(adjacencyMatrix,"block", c(.1, .2, .3, .5, .7), mySBM$memberships)
Z0 <- matrix(0, n, Q); Z0[cbind(1:n, mySBM$memberships)] <- 1
fittedSampling_block <- blockSampling_fit$new(sampledNet_block, Z0)
fittedSampling_block$vExpec
fittedSampling_block$parameters
fittedSampling_block$df
fittedSampling_block$penalty
-2 * fittedSampling_block$logLik + fittedSampling_block$penalty

sampledNet_degree <- samplingSBM(adjacencyMatrix,"degree", c(-.5,0.01))
Z0 <- matrix(0, n, Q); Z0[cbind(1:n, mySBM$memberships)] <- 1
fittedSampling_degree <- degreeSampling_fit$new(sampledNet_degree, Z0, mySBM$connectParam)
fittedSampling_degree$vExpec
fittedSampling_degree$parameters
fittedSampling_degree$df
fittedSampling_degree$penalty
-2 * fittedSampling_degree$vExpec + fittedSampling_degree$penalty

