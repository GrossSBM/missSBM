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
      link     = NULL,       # directed or not, depends on what we are modeling
      vICLs    = NULL        # ICL's of the model collection
  )
)

SBM_collection$set("public", "initialize",
  function(sample, vBlocks, sampling, family, link) {
    self$family   <- family
    self$link     <- link
    self$vBlocks  <- vBlocks
    self$vICLs    <- vector("numeric", length = length(vBlocks))
    self$sampling <- sampling
    self$sampledNetwork <- sampledNetwork$new(sample, link)
    self$samplingData   <- switch(sampling,
                            "doubleStandard" = sampling_doubleStandard$new(self$sampledNetwork$nNodes, c(.5,.5), link),
                            "class"          = sampling_class$new(self$sampledNetwork$nNodes, NA, link),
                            "starDegree"     = sampling_starDegree$new(self$sampledNetwork$nNodes, coefficients(glm(self$sampledNetwork$samplingVector~rowSums(sample, na.rm=TRUE), family = binomial(link = "logit"))), link),
                            "MAREdge"        = sampling_randomPairMAR$new(self$sampledNetwork$nNodes, .5, link),
                            "MARNode"        = sampling_randomNodeMAR$new(self$sampledNetwork$nNodes, rep(.5, self$sampledNetwork$nNodes), link),
                            "snowball"       = sampling_snowball$new(self$sampledNetwork$nNodes, rep(.5, self$sampledNetwork$nNodes), link))
  }
)

SBM_collection$set("public", "estimate",
                   function() {
                     temp <- ifelse(self$link, "Directed", "Undirected")
                     type <- paste0(self$family, temp)
                     require(parallel)
                     self$models <- lapply(self$vBlocks, function(i){
                      SBM <- switch(type,
                              "BernoulliUndirected" = SBM_BernoulliUndirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                              "BernoulliDirected"   = SBM_BernoulliDirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                              "PoissonUndirected"   = SBM_PoissonUndirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                              "PoissonDirected"     = SBM_PoissonDirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05))
                      SBM_VEMfit <- SBM_VEMfit$new(SBM, self$sampledNetwork, self$samplingData)
                      SBM_VEMfit$doVEM()
                      self$vICLs[i] <- SBM_VEMfit$vICL
                      return(SBM_VEMfit$SBM)
                       }
                      )
                   }
)

SBM_collection$set("public", "getBestModel",
                   function() {
                     if(length(which(self$vICLs == 0) > 0)){
                       return(self$models[[which.min(self$vICLs[-which(self$vICLs == 0)])]])
                     } else {
                       return(self$models[[which.min(self$vICLs)]])
                     }
                   }
)

### Tests :
# # SBM :
# mySBM <- SBM_BernoulliUndirected.fit$new(100, c(1/2, 1/2), matrix(c(.5, .05, .05, .5),2,2))
# 
# # # Sampled SBM :
# # mySampledSBM   <- sampling_doubleStandard$new(100, c(1/2, 1/2), FALSE)
# # Y              <- mySBM$rSBM()$adjacencyMatrix
# # sample         <- mySampledSBM$rSampling(Y)
# 
# # Sampled SBM 2 (MAR):
# mySampledSBM  <- sampling_randomPairMAR$new(100, 1/2, FALSE)
# Y             <- mySBM$rSBM()$adjacencyMatrix
# sampledNetwork <- mySampledSBM$rSampling(Y)
# 
# # VEM :
# sbm <- SBM_collection$new(sample$adjacencyMatrix, 2, "MAREdge", "Bernoulli", FALSE)
# sbm$estimate()
# sbm$getBestModel()
# sbm$vICLs









