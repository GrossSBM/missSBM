rm(list=ls())
library(missSBM)
library(aricode)

logistic <- function(x) {1/(1 + exp(-x))}
logit    <- function(x) {log(x/(1 - x))}

### A SBM model with covariates ###
N <- 200
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.15, Q) + .01                 # connectivity matrix
directed <- FALSE

beta  <- c(0, 1)
M <- length(beta)
covariates <- t(rmultinom(N, 1, c(1/2,1/2)))
covariates[covariates == 0] <- -1
covariates <- covariates + rnorm(N*2,0,0.1)

### Draw a undirected SBM model with covariates
mySBM <- simulateSBM(N, alpha, pi, directed, covariates, beta)

adjacencyMatrix <- mySBM$adjacencyMatrix             # the adjacency matrix

## ______________________________________________________________________
## DYAD SAMPLING

## Draw random missing entries: MAR case (dyads)
psi <- runif(ncol(covariates), -5, 5)
sampledNet <- samplingSBM(adjacencyMatrix, "dyad", psi, covariates)

cl0 <- init_clustering(sampledNet$adjacencyMatrix, Q, sampledNet$covariatesArray, "spectral")
NID(cl0,  mySBM$memberships)

## Perform inference

## control parameter for the VEM
control <- list(threshold = 1e-3, maxIter = 200, fixPointIter = 3, trace = TRUE)


missSBM <- missingSBM_fit$new(sampledNet, Q, "dyad", clusterInit = "spectral")
out_nocov <- missSBM$doVEM(control)

## Perform inference
sampledNet$covariatesArray <- covariates
missSBM_covar <- missingSBM_fit$new(sampledNet, Q, "dyad", "spectral")
out_cov <- missSBM_covar$doVEM(control)

par(mfrow = c(2,2))
plot(out_nocov$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out_nocov$objective, type = "l", main = "Variational bound along the VEM")
plot(out_cov$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out_cov$objective, type = "l", main = "Variational bound along the VEM")

print(missSBM$plot("imputedNetwork"))
print(missSBM_covar$plot("imputedNetwork"))

NID(missSBM$fittedSBM$memberships, mySBM$memberships)
NID(missSBM_covar$fittedSBM$memberships, mySBM$memberships)

print(abs(missSBM$fittedSampling$parameters - psi))
print(sum((missSBM$fittedSBM$connectParam - pi)^2))

print(abs(missSBM_covar$fittedSampling$parameters - psi))
print(sum((missSBM_covar$fittedSBM$connectParam - pi)^2))

# ## ______________________________________________________________________
# ## NODE SAMPLING
#
# ## Draw random missing entries: MAR case (node)
# psi <- 0.2
# sampledNet <- samplingSBM(adjacencyMatrix, "node", psi)
# ## Perform inference
# missSBM <- missingSBM_fit$new(sampledNet, Q, "node")
# out <- missSBM$doVEM(control)
#
# par(mfrow = c(1,2))
# plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
# plot(out$objective, type = "l", main = "Variational bound along the VEM")
#
# missSBM$plot("imputedNetwork")
#
# NID(missSBM$fittedSBM$memberships, mySBM$memberships)
# print(abs(missSBM$fittedSampling$parameters - psi))
# print(sum((missSBM$fittedSBM$connectParam - pi)^2))
