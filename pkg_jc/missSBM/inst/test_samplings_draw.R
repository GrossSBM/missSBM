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
dyad <- samplingSBM(adjacencyMatrix, "dyad", 0.1)
dyad$plot("Dyad sampling")

node <- samplingSBM(adjacencyMatrix, "node", 0.1)
node$plot("Node sampling")

block <- samplingSBM(adjacencyMatrix, "block", c(.1, .2, .7), mySBM$memberships)
block$plot("Block sampling")

double_standard <- samplingSBM(adjacencyMatrix,"double_standard", c(0.1, 0.5))
double_standard$plot("Double standard sampling")

degree <- samplingSBM(adjacencyMatrix,"degree", c(0.01,0.01))
degree$plot("Degree standard sampling")

snowball <- samplingSBM(adjacencyMatrix,"snowball", .3)
snowball$plot("Snowball standard sampling")

