library(missSBM)
library(aricode)
library(purrr)
library(sbm)
library(profvis)

## SBM parameters
N <- 10000 # number of nodes
Q <- 10   # number of clusters
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.25,Q) + .1 ) # connectivity matrix

## generate a undirected binary SBM with no covariate
sbm <- sbm::sampleSimpleSBM(N, pi, theta)
A <- sbm$networkData
diag(A) <- 0
A <- as(A, 'dgCMatrix')

vBlocks <- 1:20
clusterings <- spectral_clustering_sparse(A, vBlocks)

ARI <- clusterings %>% map( as.vector) %>% map_dbl(aricode::ARI, sbm$memberships)
plot(vBlocks, ARI, type = 'l', col = 'red')

## SBM parameters
N <- 2000 # number of nodes
Q <- 10   # number of clusters
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.25,Q) + .1 ) # connectivity matrix

## generate a undirected binary SBM with no covariate
sbm <- sbm::sampleSimpleSBM(N, pi, theta)
A <- sbm$networkData
diag(A) <- 0

vBlocks <- 1:20
clusterings <- spectral_clustering_dense(A, vBlocks)

ARI <- clusterings %>% map( as.vector) %>% map_dbl(aricode::ARI, sbm$memberships)
plot(vBlocks, ARI, type = 'l', col = 'red')
