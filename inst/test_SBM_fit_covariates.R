rm(list=ls())
library(missSBM)
library(aricode)

set.seed(1111)

### A SBM model : ###
N <- 100
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.15, Q) +.01                 # connectivity matrix
directed <- FALSE

### Draw a undirected SBM model with covariates
covarParam  <- c(-1, 0, 1)
M <- length(covarParam)
covariates <- matrix(rnorm(N*M,mean = 0, sd = 1), N, M)
mySBM <- simulateSBM(N, alpha, pi, directed, covariates, covarParam)

## testing the different initializations
## random
cat("\n VEM randomly initialized")
mySBM_fit <- SBM_fit_covariates$new(mySBM$adjacencyMatrix, mySBM$covariates, Q, sample(mySBM$memberships))
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

print(NID(mySBM_fit$memberships, mySBM$memberships))
mySBM_fit$vICL(mySBM$adjacencyMatrix)
mySBM_fit$vBound(mySBM$adjacencyMatrix)

## spectral clustering
cat("\n VEM initialized with spectral clustering")
mySBM_fit <- SBM_fit_nocovariate$new(mySBM$adjacencyMatrix, Q, "spectral")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))

## Hierarchical clustering
cat("\n VEM initialized with hierarchical clustering")
mySBM_fit <- SBM_fit_nocovariate$new(mySBM$adjacencyMatrix, Q, "hierarchical")
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)
cat("\n NID:")
print(NID(mySBM_fit$memberships, mySBM$memberships))

## Testing model selection criterion
cat("\n Assessing model selection - VEM on varying number of blocks.")
vBlocks <- 1:10
cat("\n Number of blocks =")
models <- lapply(vBlocks, function(nBlocks) {
  cat("", nBlocks)
  myFit <- SBM_fit_nocovariate$new(mySBM$adjacencyMatrix, nBlocks)
  myFit$doVEM(mySBM$adjacencyMatrix)
  myFit
})

vBound   <- sapply(models, function(model) model$vBound(mySBM$adjacencyMatrix))
vExpec   <- sapply(models, function(model) model$vExpec(mySBM$adjacencyMatrix))
vEntropy <- sapply(models, function(model) model$entropy)
vICLs    <- sapply(models, function(model) model$vICL(mySBM$adjacencyMatrix))
vPens    <- sapply(models, function(model) model$penalty)

par(mfrow=c(1,1))
plot (vBlocks,                vICLs, col = "red", type = "l")
lines(vBlocks, -2*(vBound - vEntropy), col = "green")
lines(vBlocks, -2*vBound           , col = "blue")
legend("topright",
  legend = c("vICLs = -2 * (vBound-vEntropy) + vPens", "-2* (vBound  - vEntropy)", "-2*vBound"), col=c("red", "green", "blue"),
  lty = 1, bty = "n")
abline(v = which.min(vICLs), lty="dashed")

bestICL <- models[[which.min(vICLs)]]

