rm(list=ls())
library(missSBM)
library(aricode)
set.seed(78910)

### A SBM model : ###
N <- 200
Q <- 3
alpha <- rep(1,Q)/Q                     # mixture parameter
pi <- diag(.15, Q) + .01                 # connectivity matrix
directed <- FALSE

covarParam  <- c(0, 1)
M <- length(covarParam)
X <- t(rmultinom(N, 1, c(1/2,1/2)))
X[X == 0] <- -1
X <- X + rnorm(N*2,0,0.1)
covariates <- array(dim = c(N, N, M))
for (i in 1:N)
  for (j in 1:N)
    covariates[i,j,] <- -abs(X[i, ] - X[j, ])

logistic <- function(x) {1/(1 + exp(-x))}
logit    <- function(x) {log(x/(1 - x))}

hist(logistic(X %*% covarParam))

mySBM <- simulateSBM(N, alpha, pi, directed, covariates, covarParam)
plot(mySBM)

## testing the different initializations
## spectral clustering
cat("\n VEM initialized with hierarchical clustering on residuals")
clusterInit <- init_clustering(mySBM$adjacencyMatrix, Q, covariates, "spectral")
mySBM_fit <- SBM_fit_covariates$new(mySBM$adjacencyMatrix, mySBM$covariates, clusterInit)
out <- mySBM_fit$doVEM(mySBM$adjacencyMatrix, trace = TRUE)

par(mfrow = c(1,2))
plot(out$delta    , type = "l", main = "Variations of SBM parameters along the VEM")
plot(out$objective, type = "l", main = "Variational bound along the VEM")

print(NID(mySBM_fit$memberships, mySBM$memberships))
mySBM_fit$vICL(mySBM$adjacencyMatrix)
mySBM_fit$vBound(mySBM$adjacencyMatrix)

## Testing model selection criterion
cat("\n Assessing model selection - VEM on varying number of blocks.")
vBlocks <- 1:5
cat("\n Number of blocks =")
models <- lapply(vBlocks, function(nBlocks) {
  cat("", nBlocks)
  clusterInit <- init_clustering(mySBM$adjacencyMatrix, nBlocks, covariates, "spectral")
  myFit <- SBM_fit_covariates$new(mySBM$adjacencyMatrix, mySBM$covariates, clusterInit)
  myFit$doVEM(mySBM$adjacencyMatrix)
  myFit
})

vBound   <- sapply(models, function(model) model$vBound(mySBM$adjacencyMatrix))
vExpec   <- sapply(models, function(model) model$vExpec(mySBM$adjacencyMatrix))
vEntropy <- sapply(models, function(model) model$entropy)
vICLs    <- sapply(models, function(model) model$vICL(mySBM$adjacencyMatrix))
vPens    <- sapply(models, function(model) model$penalty)

par(mfrow=c(1,1))
plot (vBlocks,                vICLs, col = "red", type = "l", ylim = range(c(vICLs, -2*vBound)))
lines(vBlocks, -2*(vBound - vEntropy), col = "green")
lines(vBlocks, -2*vBound           , col = "blue")
legend("topright",
  legend = c("vICLs = -2 * (vBound-vEntropy) + vPens", "-2* (vBound  - vEntropy)", "-2*vBound"), col=c("red", "green", "blue"),
  lty = 1, bty = "n")
abline(v = which.min(vICLs), lty = "dashed")

bestICL <- models[[which.min(vICLs)]]

logistic(bestICL$connectParam)

