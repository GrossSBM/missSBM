library(R6)
source("R/utils.R")
source("R/SBM-R6Class.R")
source("R/SBMfit-R6Class.R")
## test between old and new SBM
n <- 100
Q <- 3
alpha <- rep(1,Q)/Q                                                                # mixture parameter
pi <- diag(.45,Q) + .05                                                            # connectivity matrix
family <- "Bernoulli"                                                              # the emission law
directed <- FALSE
## tested for Bernoulli, Poisson, directed or not : works

mySBM_old <- SBM_BernoulliUndirected$new(n, alpha, pi)
mySBM_old$rBlocks()
mySBM_old$rAdjMatrix()

mySBM_new <- SBM_fit$new(mySBM_old$adjacencyMatrix, Q, clusterInit = "hierarchical")
print(mySBM_new$cLogLik())

mySBM_new$blockVarPar <- mySBM_new$blocks
mySBM_new$fixPoint()

# tau <- mySBM_new$blockVarPar
# X <- mySBM_new$adjacencyMatrix
# Z <- mySBM_new$blocks
# edges <- matrix(TRUE, n, n); diag(edges) <- FALSE
# if (!directed) edges[lower.tri(edges)] <- FALSE
#
# d_law <- function(x, prob) {dbinom(x, 1, prob)}
#
# term1 <- X %*% tau %*% t(log(pi)) + bar(X) %*% tau %*% t(log(1 - pi))

# term2 <- t(log(d_law(X, tau %*% pi %*% t(tau)) %*% tau )

# sum( log (
#   d_law(outer(X[1, ], rep(1,Q)), outer(rep(1,n), pi[1,]))
#   ) * c(0,rep(1,n-1)) * tau
# )
#
# In <- matrix(1, n, n)
# diag(In) <- 0
# mask <- rep(1,Q) %x% In
#
# colSums( as.vector(tau) * (log (
#   d_law(kronecker(rep(1,Q), X), kronecker(pi[1, ], rep(1,n)))
#   ) * mask )
# )
#
# colSums( as.vector(tau) * (log (
#   d_law(kronecker(rep(1,Q), X), kronecker(pi[1,], rep(1,n)))
#   ) * mask )
# )
