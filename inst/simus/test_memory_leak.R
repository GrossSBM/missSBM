library(missSBM)
library(sbm)

## SBM parameters
N <- 400 # number of nodes
Q <- 5   # number of clusters
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.25,Q) + .1 ) # connectivity matrix

## generate a undirected binary SBM with no covariate
sbm <- sbm::sampleSimpleSBM(N, pi, theta)
A <- sbm$networkData
diag(A) <- 0

vBlocks <- 1:10
clusterings <- spectral_clustering_dense(A, vBlocks)
