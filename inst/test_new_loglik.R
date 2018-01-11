library(R6)
source("R/SBM-R6Class.R")
source("R/SBMfit-R6Class.R")
## test between old and new SBM
n <- 300
Q <- 3
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emission law

## tested for Bernoulli, Poisson, directed or not : works

mySBM_old <- SBM_PoissonUndirected$new(n, alpha, pi)
mySBM_old$rBlocks()
mySBM_old$rAdjMatrix()

mySBM_new <- SBM_fit$new("Poisson", FALSE, n, alpha, pi)
mySBM_new$blocks <- mySBM_old$blocks
mySBM_new$adjacencyMatrix <- mySBM_old$adjacencyMatrix

print(mySBM_old$cLogLik)

print(mySBM_new$cLogLik())
