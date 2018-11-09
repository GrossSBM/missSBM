rm(list=ls())
library(missSBM)
library(aricode)

### A SBM model : ###
n <- 400
Q <- 5
alpha <- rep(1,Q)/Q       # mixture parameter
pi <- diag(.45,Q) + .05   # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(n, alpha, pi, directed) # simulation of a Bernoulli non-directed SBM
adjacencyMatrix <- mySBM$adjacencyMatrix             # the adjacency matrix

## Draw random missing entries: MAR case (dyads)
psi <- 0.3
sampledNet <- samplingSBM(adjacencyMatrix, "dyad", psi)
A_dyad <- sampledNet$adjacencyMatrix

psi <- 0.3
sampledNet <- samplingSBM(adjacencyMatrix, "node", psi)
A_node <- sampledNet$adjacencyMatrix

NID(init_clustering(A_dyad, Q, NULL, "spectral"    ), mySBM$memberships)
NID(init_clustering(A_dyad, Q, NULL, "kmeans"      ), mySBM$memberships)
NID(init_clustering(A_dyad, Q, NULL, "hierarchical"), mySBM$memberships)

NID(init_clustering(A_node, Q, NULL, "spectral"    ), mySBM$memberships)
NID(init_clustering(A_node, Q, NULL, "kmeans"      ), mySBM$memberships)
NID(init_clustering(A_node, Q, NULL, "hierarchical"), mySBM$memberships)

