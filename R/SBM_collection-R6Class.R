#' a collection of adjusted Stochastic Block Model
#'
#' @importFrom R6 R6Class
SBM_collection <-
  R6Class(classname = "SBM_collection",
    public = list(
      sampledNetwork = NULL, # the sampled network data
      vBlocks  = NULL,       # the vector of number of blocks considered in the collection
      models   = NULL,       # the collection of SBM fit
      sampling = NULL,       # the sampling design for missing data modeling
      samplingData = NULL,   # informations about the sampling
      family   = NULL,       # the emission law of the adjacency matrix
      link     = NULL        # directed or not, depends on what we are modeling
  )
)

SBM_collection$set("public", "initialize",
  function(sampledNetwork, vBlocks, sampling, family, link) {
    self$family   <- family
    self$link     <- link
    self$vBlocks  <- vBlocks
    self$sampling <- sampling
    self$sampledNetwork <- sampledNetwork$new(sampledNetwork, link)
    self$samplingData   <- switch(sampling,
                            "doubleStandard" = sampling_doubleStandard$new(self$sampledNetwork$nNodes, NA, link),
                            "class"          = sampling_class$new(self$sampledNetwork$nNodes, NA, link),
                            "starDegree"     = sampling_starDegree$new(self$sampledNetwork$nNodes, NA, link),
                            "MAREdge"        = sampling_randomPairMAR$new(self$sampledNetwork$nNodes, NA, link),
                            "MARNode"        = sampling_randomNodeMAR$new(self$sampledNetwork$nNodes, NA, link),
                            "snowball"       = sampling_snowball$new(self$sampledNetwork$nNodes, NA, link))
  }
)

SBM_collection$set("public", "estimate",
                   function() {
                     temp <- ifelse(link, "Directed", "Undirected")
                     type <- paste0(family, temp)
                     require(parallel)
                     self$models <- mclapply(vBlocks, function(i){
                      SBM <- switch(type,
                              "BernoulliUndirected" = SBM_BernoulliUndirected.fit$new(self$sampledNetwork$nNodes, rep(0, i), NA),
                              "BernoulliDirected"   = SBM_BernoulliDirected.fit$new(self$sampledNetwork$nNodes, rep(0, i), NA),
                              "PoissonUndirected"   = SBM_PoissonUndirected.fit$new(self$sampledNetwork$nNodes, rep(0, i), NA),
                              "PoissonDirected"     = SBM_PoissonDirected.fit$new(self$sampledNetwork$nNodes, rep(0, i), NA))
                      SBM_VEMfit <- SBM_VEMfit$new(SBM, self$sampledNetwork, self$samplingData)
                      SBM_VEMfit$doVEM()
                      return(SBM)
                       }
                     ,mc.cores = detectCores()-1)
                     SBM_VEMfit$new(mySBM, sampledNetwork, mySampledSBM)
                   }
)

### Tests :

# SBM : 
mySBM <- SBM_BernoulliUndirected.fit$new(100, c(1/2, 1/2), matrix(c(.5, .05, .05, .5),2,2))

# Sampled SBM :
mySampledSBM   <- sampling_doubleStandard$new(100, c(1/2, 1/2), FALSE)
Y              <- mySBM$rSBM()$adjacencyMatrix
sampledNetwork <- mySampledSBM$rSampling(Y)

# # Sampled SBM 2 (MAR):
# mySampledSBM  <- sampling_randomPairMAR$new(100, 1/2, FALSE)
# Y             <- mySBM$rSBM()$adjacencyMatrix
# sampledNetwork <- mySampledSBM$rSampling(Y)

# VEM : 
sbm <- SBM_collection$new(sampledNetwork, 1:5, "doubleStandard", "Bernoulli", FALSE)











