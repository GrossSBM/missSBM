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

## contorl parameter for the VEM
control <- list(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = TRUE)

## ______________________________________________________________________
## DYAD SAMPLING

## Draw random missing entries: MAR case (dyads)
psi <- 0.3
sampledNet <- samplingSBM(adjacencyMatrix, "dyad", psi)
## Perform inference
missSBM <- missingSBM_fit$new(sampledNet, Q, "dyad")
out <- missSBM$doVEM(control)


par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

print(missSBM$plot("imputedNetwork"))

NID(missSBM$fittedSBM$memberships, mySBM$memberships)
print(abs(missSBM$fittedSampling$parameters - psi))
print(sum((missSBM$fittedSBM$connectParam - pi)^2))

## ______________________________________________________________________
## NODE SAMPLING

## Draw random missing entries: MAR case (node)
psi <- 0.2
sampledNet <- samplingSBM(adjacencyMatrix, "node", psi)
## Perform inference
missSBM <- missingSBM_fit$new(sampledNet, Q, "node")
out <- missSBM$doVEM(control)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

missSBM$plot("imputedNetwork")

NID(missSBM$fittedSBM$memberships, mySBM$memberships)
print(abs(missSBM$fittedSampling$parameters - psi))
print(sum((missSBM$fittedSBM$connectParam - pi)^2))

## ______________________________________________________________________
## DOUBLE-STANDARD SAMPLING

## Draw random missing entries: NMAR case (double_standard)
psi <- c(.3, .6)
sampledNet <- samplingSBM(adjacencyMatrix, "double_standard", psi)
## Perform inference
missSBM <- missingSBM_fit$new(sampledNet, Q, "double_standard")
out <- missSBM$doVEM(control)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

print(missSBM$plot("imputedNetwork"))

NID(missSBM$fittedSBM$memberships, mySBM$memberships)
print(abs(missSBM$fittedSampling$parameters - psi))
print(sum((missSBM$fittedSBM$connectParam - pi)^2))

## ______________________________________________________________________
## BLOCK SAMPLING

## Draw random missing entries: NMAR case (blocks)
psi <- c(.1, .3, .2, .5, .7)
sampledNet <- samplingSBM(adjacencyMatrix, "block", psi, mySBM$memberships)
## Perform inference
missSBM <- missingSBM_fit$new(sampledNet, Q, "block")
out <- missSBM$doVEM(control)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

print(missSBM$plot("imputedNetwork"))

NID(missSBM$fittedSBM$memberships, mySBM$memberships)
print(abs(sort(missSBM$fittedSampling$parameters) - sort(psi)))
print(sum((missSBM$fittedSBM$connectParam - pi)^2))

## ______________________________________________________________________
## DEGREE SAMPLING

## Draw random missing entries: NMAR case (blocks)
psi <- c(-5, .1)
sampledNet <- samplingSBM(adjacencyMatrix, "degree", psi)
## Perform inference
missSBM <- missingSBM_fit$new(sampledNet, 5, "degree")
out <- missSBM$doVEM(control)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

missSBM$plot("imputedNetwork")

NID(missSBM$fittedSBM$memberships, mySBM$memberships)
print(abs(missSBM$fittedSampling$parameters - psi))
print(sum((missSBM$fittedSBM$connectParam - pi)^2))
missSBM$fittedSBM$connectParam

logistic <- function(x) {1/(1 + exp(-x))}

plot(logistic(psi[1] + psi[2] * rowSums(adjacencyMatrix)),
     logistic(missSBM$fittedSampling$parameters[1] + missSBM$fittedSampling$parameters[2] * rowSums(missSBM$imputedNetwork)))
abline(0,1)
