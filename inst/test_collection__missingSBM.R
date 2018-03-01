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
mySBM <- simulateSBM(n, alpha, pi, directed) # simulation of ad Bernoulli non-directed SBM
adjacencyMatrix <- mySBM$adjacencyMatrix             # the adjacency matrix

## ______________________________________________________________________
## DYAD SAMPLING

## Draw random missing entries: MAR case (dyads)
psi <- 0.3
sampledNet <- samplingSBM(adjacencyMatrix, "dyad", psi)

vBlocks <- 1:8
out <- inferSBM(sampledNet$adjacencyMatrix, vBlocks, "dyad")
models <- out$models
vpens <- sapply(models, function(model) model$penalty)
ventr <- sapply(models, function(model) model$entropyImputed)
vBICs <- sapply(models, function(model) model$vBIC)
vICLs <- sapply(models, function(model) model$vICL)
vJ <- sapply(out$models, function(model) model$vBound)
par(mfrow=c(2,2))
plot(vBlocks, vICLs, type = "l")
plot(vBlocks, vJ, type = "l")
plot(vBlocks, ventr, type = "l")
plot(vBlocks, vpens, type = "l")
#plot(out$monitor$objective[out$monitor$nBlocks == Q], type = "l")

best <- models[[which.min(vICLs)]]
best$plot()

## ______________________________________________________________________
## NODE SAMPLING

## Draw random missing entries: MAR case (nodes)
psi <- 0.3
sampledNet <- samplingSBM(adjacencyMatrix, "node", psi)

vBlocks <- 1:8
out <- inferSBM(sampledNet$adjacencyMatrix, vBlocks, "node")
models <- out$models
vJ_SBM <- sapply(models, function(model) model$fittedSBM$vBound(model$imputedNetwork))
vJ_sampling <- sapply(models, function(model) model$fittedSampling$logLik)
vpens <- sapply(models, function(model) model$penalty)
ventr <- sapply(models, function(model) model$entropyImputed)
vBICs <- sapply(models, function(model) model$vBIC)
vICLs <- sapply(models, function(model) model$vICL)
vJ <- sapply(out$models, function(model) model$vBound)
par(mfrow=c(2,2))
plot(vBlocks, vICLs, type = "l")
plot(vBlocks, vJ, type = "l")
plot(vBlocks, ventr, type = "l")
plot(vBlocks, vpens, type = "l")
best <- models[[which.min(vICLs)]]
best$plot()

## ______________________________________________________________________
## DOUBLE-STANDARD SAMPLING

## Draw random missing entries: NMAR case (double_standard)
psi <- c(.3, .6)
sampledNet <- samplingSBM(adjacencyMatrix, "double_standard", psi)

vBlocks <- 1:8
out <- inferSBM(sampledNet$adjacencyMatrix, vBlocks, "double_standard")
models <- out$models
vICLs <- sapply(models, function(model) model$vICL)
vJ <- sapply(out$models, function(model) model$vBound)
par(mfrow=c(1,2))
plot(vBlocks, vICLs, type = "l")
plot(vBlocks, vJ, type = "l")
best <- models[[which.min(vICLs)]]
best$plot()

## ______________________________________________________________________
## BLOCK SAMPLING

## Draw random missing entries: NMAR case (double_standard)
psi <- c(.1, .3, .2, .5, .7)
sampledNet <- samplingSBM(adjacencyMatrix, "block", psi, mySBM$memberships)

vBlocks <- 1:10
out <- inferSBM(sampledNet$adjacencyMatrix, vBlocks, "block")
models <- out$models
vICLs <- sapply(models, function(model) model$vICL)
vJ <- sapply(out$models, function(model) model$vBound)
par(mfrow=c(1,2))
plot(vBlocks, vICLs, type = "l")
plot(vBlocks, vJ, type = "l")
best <- models[[which.min(vICLs)]]
best$plot()

