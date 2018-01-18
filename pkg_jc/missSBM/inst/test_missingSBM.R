library(missSBM)
library(aricode)

### A SBM model : ###
n <- 300
Q <- 3
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emmission law
directed <- FALSE                                                                  # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM
adjacencyMatrix <- mySBM$adjacencyMatrix                                           # the adjacency matrix

## Draw random missing entries: MAR case (dyads)
psi <- 0.1
sampledNet <- samplingSBM(adjacencyMatrix, "dyad", psi)
## Perform inference
obj <- missingSBM_fit$new(sampledNet$adjacencyMatrix, 3, "dyad")
out <- obj$doVEM(trace = TRUE)
NID(obj$fittedSBM$memberships, mySBM$memberships)
print(abs(obj$fittedSampling$parameters - psi))
print(sum((obj$fittedSBM$connectParam - pi)^2))

## Draw random missing entries: MAR case (node)
psi <- 0.1
sampledNet <- samplingSBM(adjacencyMatrix, "node", psi)
## Perform inference
obj <- missingSBM_fit$new(sampledNet$adjacencyMatrix, 3, "node")
out <- obj$doVEM(trace = TRUE)
NID(obj$fittedSBM$memberships, mySBM$memberships)
print(abs(obj$fittedSampling$parameters - psi))
print(sum((obj$fittedSBM$connectParam - pi)^2))

## Draw random missing entries: NMAR case (blocks)
psi <- c(.1, .2, .7)
sampledNet <- samplingSBM(adjacencyMatrix, "block", psi, mySBM$memberships)
## Perform inference
obj <- missingSBM_fit$new(sampledNet$adjacencyMatrix, 3, "block")
out <- obj$doVEM(trace = TRUE)
NID(obj$fittedSBM$memberships, mySBM$memberships)
print(abs(obj$fittedSampling$parameters - psi))
print(sum((obj$fittedSBM$connectParam - pi)^2))

## Draw random missing entries: NMAR case (double_standard)
psi <- c(.3, .6)
sampledNet <- samplingSBM(adjacencyMatrix, "double_standard", psi)
## Perform inference
obj <- missingSBM_fit$new(sampledNet$adjacencyMatrix, 3, "double_standard")
out <- obj$doVEM(trace = TRUE)
NID(obj$fittedSBM$memberships, mySBM$memberships)
print(abs(obj$fittedSampling$parameters - psi))
print(sum((obj$fittedSBM$connectParam - pi)^2))


