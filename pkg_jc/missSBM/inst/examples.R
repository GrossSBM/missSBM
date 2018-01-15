### A SBM model : ###
n <- 300
Q <- 3
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emmission law
directed <- FALSE                                                                  # if the network is directed or not
mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM

### Results : ###
adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix

## Sampling of the data : ##
samplingParameters <- .5                                                           # the sampling rate
sampling <- "MAREdge"                                                              # the sampling design
sampledAdjMatrix <- samplingSBM(adjacencyMatrix, sampling, samplingParameters)     # the sampled adjacency matrix


# ## Inference :
# vBlocks <- 1:5                                                                     # number of classes
# sbm <- inferSBM(sampledAdjMatrix, vBlocks, sampling)             # the inference
