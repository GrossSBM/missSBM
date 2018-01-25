library(missSBM)
library(aricode)

set.seed(1111)

### A SBM model : ###
n <- 300
Q <- 3
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emmission law
directed <- FALSE                                                                  # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM

## testing the different initializations
## random
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, 3, sample(mySBM$memberships))
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

print(NID(mySBM_fit$memberships, mySBM$memberships))
mySBM_fit$vICL(mySBM$adjacencyMatrix)
mySBM_fit$vLogLik(mySBM$adjacencyMatrix)
mySBM_fit$vBIC(mySBM$adjacencyMatrix)

## spectral clustering
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, 3, "spectral")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))

## Hierarchical clustering
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, 3, "hierarchical")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))

## K-means clustering
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, 3, "kmeans")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))


