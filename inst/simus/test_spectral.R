library(missSBM)
library(aricode)
library(purrr)
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
A <- as(A, 'dgCMatrix')
myNet <- missSBM:::partlyObservedNetwork$new(A)

vBlocks <- 1:10
timings_new <- system.time(clusterings_new <- myNet$clustering_new(vBlocks))[3]
timings_old <- system.time(clusterings_old <- myNet$clustering(vBlocks))[3]

ARI_new <- clusterings_new %>% map( as.vector) %>% map_dbl(ARI, sbm$memberships)
ARI_old <- clusterings_old %>% map( as.vector) %>% map_dbl(ARI, sbm$memberships)

plot(vBlocks, ARI_old, type = 'l', col = 'red')
lines(vBlocks, ARI_new, type = 'l', col = 'blue')
