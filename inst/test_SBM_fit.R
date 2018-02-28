library(missSBM)
library(aricode)

set.seed(1111)

### A SBM model : ###
n <- 400
Q <- 5
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emmission law
directed <- FALSE                                                                  # if the network is directed or not

### Draw a SBM model
mySBM <- simulateSBM(n, alpha, pi, family, directed)                               # simulation of ad Bernoulli non-directed SBM

## testing the different initializations
## random
cat("\n VEM randomly initialized")
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, Q, sample(mySBM$memberships))
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

print(NID(mySBM_fit$memberships, mySBM$memberships))
mySBM_fit$vICL(mySBM$adjacencyMatrix)
mySBM_fit$vBound(mySBM$adjacencyMatrix)
mySBM_fit$vBIC(mySBM$adjacencyMatrix)

## spectral clustering
cat("\n VEM initialized with spectral clustering")
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, Q, "spectral")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))

## Hierarchical clustering
cat("\n VEM initialized with hierarchical clustering")
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, Q, "hierarchical")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))

## K-means clustering
cat("\n VEM initialized with K-means clustering")
mySBM_fit <- SBM_fit$new(mySBM$adjacencyMatrix, Q, "kmeans")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))

## Testing model selection criterion
cat("\n Assessing model selection - VEM on varying number of blocks.")
vBlocks <- 1:10
cat("\n Number of blocks =")
models <- lapply(vBlocks, function(nBlocks) {
  cat("", nBlocks)
  myFit <- SBM_fit$new(mySBM$adjacencyMatrix, nBlocks)
  myFit$doVEM(mySBM$adjacencyMatrix)
  myFit
})

vJs   <- sapply(models, function(model) model$vBound(mySBM$adjacencyMatrix))
vICLs <- sapply(models, function(model) model$vICL(mySBM$adjacencyMatrix))
vBICs <- sapply(models, function(model) model$vBIC(mySBM$adjacencyMatrix))
par(mfrow = c(1,3))
plot(vBlocks, -2*vJs, type = "l")
plot(vBlocks, vBICs , type = "l")
plot(vBlocks, vICLs , type = "l")
bestICL <- models[[which.min(vICLs)]]

