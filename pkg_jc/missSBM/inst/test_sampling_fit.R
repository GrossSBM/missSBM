library(missSBM)

### A SBM model : ###
n <- 300
Q <- 3
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emmission law
directed <- FALSE                                                                  # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM

adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix

## Draw
sampledNet_dyad <- samplingSBM(adjacencyMatrix, "dyad", 0.1)
fittedSampling_dyad <- dyadSampling_fit$new(sampledNet_dyad$adjacencyMatrix)
fittedSampling_dyad$vLogLik
fittedSampling_dyad$parameters
fittedSampling_dyad$df
fittedSampling_dyad$penalty
fittedSampling_dyad$vLogLik + fittedSampling_dyad$penalty

sampledNet_node <- samplingSBM(adjacencyMatrix, "node", 0.1)
fittedSampling_node <- nodeSampling_fit$new(sampledNet_node$adjacencyMatrix)
fittedSampling_node$vLogLik
fittedSampling_node$parameters
fittedSampling_node$df
fittedSampling_node$penalty
fittedSampling_node$vLogLik + fittedSampling_dyad$penalty

sampledNet_double_standard <- samplingSBM(adjacencyMatrix,"double_standard", c(0.1, 0.5))
fittedSampling_double_standard <- doubleStandardSampling_fit$new(sampledNet_double_standard$adjacencyMatrix)
fittedSampling_double_standard$vLogLik
fittedSampling_double_standard$parameters
fittedSampling_double_standard$df
fittedSampling_double_standard$penalty
fittedSampling_double_standard$vLogLik + fittedSampling_dyad$penalty

sampledNet_block <- samplingSBM(adjacencyMatrix,"block", c(.1, .2, .7), mySBM$memberships)
Z0 <- matrix(0, n, Q); Z0[cbind(1:n, mySBM$memberships)] <- 1
fittedSampling_block <- blockSampling_fit$new(sampledNet_block$adjacencyMatrix, Z0)
fittedSampling_block$vLogLik
fittedSampling_block$parameters
fittedSampling_block$df
fittedSampling_block$penalty
fittedSampling_block$vLogLik + fittedSampling_dyad$penalty


sampledNet_degree <- samplingSBM(adjacencyMatrix,"degree", c(0.01,0.01))
fittedSampling_degree <- degreeSampling_fit$new(sampledNet_degree$adjacencyMatrix)
fittedSampling_degree$vLogLik
fittedSampling_degree$parameters
fittedSampling_degree$df
fittedSampling_degree$penalty
fittedSampling_degree$vLogLik + fittedSampling_dyad$penalty

#
# snowball <- samplingSBM(adjacencyMatrix,"snowball", .3)

