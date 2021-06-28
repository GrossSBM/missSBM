library(missSBM)
library(aricode)
library(sbm)
library(Matrix)

## SBM parameters
N <- 2000 # number of nodes
Q <- 5   # number of clusters
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.25,Q) + .1 ) # connectivity matrix

## generate a undirected binary SBM with no covariate
sbm <- sbm::sampleSimpleSBM(N, pi, theta)
A <- sbm$networkData
diag(A) <- 0
A <- Matrix(A)
out <- missSBM:::spectral_clustering(A, 10)

Un <- out[, 1:5]
Un <- sweep(Un, 1, sqrt(rowSums(Un^2)), "/")
pairs(Un)

ARI(kmeans(Un, 5)$cl, sbm$memberships)

