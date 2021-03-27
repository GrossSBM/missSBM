library(sbm)
library(missSBM)

set.seed(178303)
### A SBM model : ###
N <- 1000
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- list(mean = diag(.45, Q, Q) + .05) # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a undirected SBM model
mySBM <- sbm::sampleSimpleSBM(N, pi, theta)
A <- observeNetwork(mySBM$networkData, "dyad", parameters = 0.5)
obsNet <- missSBM:::partlyObservedNetwork$new(A)
cl <- obsNet$clustering(nbBlocks = Q, "spectral")

res <- microbenchmark::microbenchmark(
  MAR = {
    mySBM_fit <- missSBM:::SimpleSBM_fit_MAR$new(A, cl)
    optim <- mySBM_fit$doVEM(trace = FALSE)
  },
  missSBM = {
    my_missSBM_fit <- missSBM:::missSBM_fit$new( missSBM:::partlyObservedNetwork$new(A), Q, "dyad", cl, FALSE)
    my_missSBM_fit$doVEM(control = list(threshold = 1e-3, maxIter = 100, fixPointIter = 5, trace = 0))
  },
  times = 10
)

ggplot2::autoplot(res)

my_missSBM_fit$fittedSBM$connectParam
mySBM_fit$connectParam
mySBM$connectParam
