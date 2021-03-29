library(sbm)
library(missSBM)

## Common parameters
nbNodes  <- 500
nbBlocks <- 2
blockProp <- c(.5, .5) # group proportions
covarParam <- c(-2,2)
dimLabels <- list(row = "rowLabel", col = "colLabel")
covar1 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covar2 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covarList_directed <- list(covar1 = covar1, covar2 = covar2)

covar1 <- covar1 + t(covar1)
covar2 <- covar2 + t(covar2)
covarList <- list(covar1 = covar1, covar2 = covar2)

means <- diag(.4, 2) + 0.05
connectParam <- list(mean = means)

## Basic construction - check for wrong specifications
mySampler <- SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList[1])
mySampler$rMemberships(store = TRUE)
mySampler$rEdges(store = TRUE)
adjacencyMatrix <- missSBM::observeNetwork(mySampler$networkData, "dyad", 1)


Y <- adjacencyMatrix
Y[is.na(Y)] <- 0
Z <- mySampler$indMemberships
T <- matrix(0.5, nrow(Z), ncol(Z))
M <- mySampler$covarEffect
Gamma <- missSBM:::.logit(mySampler$connectParam$mean)
pi <- mySampler$blockProp

# Storing data
## where are my observations?
diag(adjacencyMatrix) <- NA
obs <- which(!is.na(adjacencyMatrix), arr.ind = TRUE)
R_sp <- Matrix::sparseMatrix(obs[,1], obs[,2],x = 1, dims = dim(adjacencyMatrix))

## where are my non-zero entries?
nzero <- which(!is.na(adjacencyMatrix) & adjacencyMatrix != 0, arr.ind = TRUE)
Y_sp   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(adjacencyMatrix))

microbenchmark::microbenchmark(
 spmat = missSBM:::vLL_complete_sparse_bernoulli_undirected_covariates(Y_sp, R_sp, M, T, Gamma, pi),
 dense = missSBM:::vExpec_covariates(Y, M, Gamma, T, pi)
) -> res

ggplot2::autoplot(res)




