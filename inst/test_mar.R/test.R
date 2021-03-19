library(aricode)

set.seed(178303)
### A SBM model : ###
N <- 200
Q <- 3
pi <- rep(1, Q)/Q           # block proportion
theta <- list(mean = diag(.45, Q, Q) + .05) # connectivity matrix
directed <- FALSE         # if the network is directed or not

### Draw a undirected SBM model
mySBM <- sbm::sampleSimpleSBM(N, pi, theta)
A <- observeNetwork(mySBM$networkData, "dyad", parameters = 0.5)
obsNet <- missSBM:::partlyObservedNetwork$new(A)
cl <- obsNet$clustering(nbBlocks = Q)

mySBM_fit <- missSBM:::SimpleSBM_fit_MAR$new(A, cl)
optim <- mySBM_fit$doVEM()
mySBM_fit$connectParam

