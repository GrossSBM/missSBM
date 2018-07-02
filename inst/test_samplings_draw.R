rm(list = ls())
library(missSBM)

### A SBM model : ###
N <- 300
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.45,Q) + .05                 # connectivity matrix
directed <- FALSE

### Draw a SBM model (Bernoulli, undirected)
mySBM <- simulateSBM(N, alpha, pi, directed)

adjacencyMatrix <- mySBM$adjacencyMatrix

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

### Draw a SBM model (Bernoulli, undirected) with covariates
M <- 10
covariates <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
covarParam  <- rnorm(M,0,1)
mySBM <- simulateSBM(N, alpha, pi, directed, covariates, covarParam)

adjacencyMatrix <- mySBM$adjacencyMatrix

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

