#' a collection of adjusted Stochastic Block Model
#'
#' @importFrom R6 R6Class
SBM_collection <-
  R6Class(classname = "SBM_collection",
    public = list(
      sampledNetwork = NULL, # the sampled network data
      vBlocks  = NULL, # the vector of number of blocks considered in the collection
      models   = NULL, # the collection of SBM fit
      sampling = NULL, # the sampling design for missing data modeling
      samplingData = NULL, # informations about the sampling
      family   = NULL, # the emission law of the adjacency matrix
      link     = NULL  # directed or not, depends on what we are modeling
      ## methods
      ### TODO: VE-step and M-step
      ## criteria = NULL, # variational bound and ICL for all models
  )
)

SBM_collection$set("public", "initialize",
  function(sampledNetwork, vBlocks, sampling, family, link) {
    self$family   <- family
    self$link     <- link
    self$vBlocks  <- vBlocks
    self$sampledNetwork <- sampledNetwork$new(sampledNetwork, link)
    self$sampling <- sampling
    self$samplingData       <- switch(sampling,
                            "doubleStandard" = sampling_doubleStandard$new(self$sampledNetwork$nNodes, NA, link),
                            "class"          = sampling_class$new(self$sampledNetwork$nNodes, NA, link),
                            "starDegree"     = sampling_starDegree$new(self$sampledNetwork$nNodes, NA, link),
                            "MAREdge"        = sampling_randomPairMAR$new(self$sampledNetwork$nNodes, NA, link),
                            "MARNode"        = sampling_randomNodeMAR$new(self$sampledNetwork$nNodes, NA, link),
                            "snowball"       = sampling_snowball$new(self$sampledNetwork$nNodes, NA, link))
  
  }
)

SBM, sampledNetwork, sampling, nBlocks