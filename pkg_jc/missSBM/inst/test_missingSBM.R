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
sampledNet_dyad <- samplingSBM(adjacencyMatrix, "dyad", 0.6)
obj <- missingSBM_fit$new(sampledNet_dyad$adjacencyMatrix, 3, "dyad")
obj$doVEM()
