library(missSBM)
library(sbm)
library(aricode)
library(profvis)

## SBM parameters
N <- 10000 # number of nodes
Q <- 5   # number of clusters
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.25,Q) + .1 ) # connectivity matrix

## generate a undirected binary SBM with no covariate
sbm <- sbm::sampleSimpleSBM(N, pi, theta)

prof_total <- profvis({
  networkData <- missSBM:::partlyObservedNetwork$new(sbm$networkData)
  cl0 <- networkData$clustering(5)[[1]]
  mySBM <- missSBM:::SimpleSBM_fit_noCov$new(networkData, cl0)
  mySBM$doVEM(threshold = 1e-5, maxIter = 50, fixPointIter = 10)
})

ARI(mySBM$memberships, sbm$memberships)

sum( (mySBM$connectParam$mean - sbm$connectParam$mean)^2)

htmlwidgets::saveWidget(prof_total, "profile.html")

prof_imputation <- profvis({
  networkData <- missSBM:::partlyObservedNetwork$new(sbm$networkData)
  A <- networkData$imputation()
})

htmlwidgets::saveWidget(prof_imputation, "~/Desktop/prof_imput.html")
