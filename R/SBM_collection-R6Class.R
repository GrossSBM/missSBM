#' A collection of adjusted Stochastic Block Model
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
      vICLs    = NULL        # ICL's of the model collection
  )
)

SBM_collection$set("public", "initialize",
  function(sample, vBlocks, sampling, family, directed) {
    self$family   <- family
    self$vBlocks  <- vBlocks
    self$sampling <- sampling
    self$sampledNetwork <- sampledNetwork$new(sample, directed)
    self$samplingData   <- switch(sampling,
                            "doubleStandard" = sampling_doubleStandard$new(self$sampledNetwork$nNodes, c(.5,.5), directed),
                            "class"          = sampling_class$new(self$sampledNetwork$nNodes, .5, directed),
                            "starDegree"     = sampling_starDegree$new(self$sampledNetwork$nNodes, coefficients(glm(self$sampledNetwork$samplingVector~rowSums(sample, na.rm=TRUE), family = binomial(directed = "logit"))), directed),
                            "MAREdge"        = sampling_randomPairMAR$new(self$sampledNetwork$nNodes, .5, directed),
                            "MARNode"        = sampling_randomNodeMAR$new(self$sampledNetwork$nNodes, rep(.5, self$sampledNetwork$nNodes), directed),
                            "snowball"       = sampling_snowball$new(self$sampledNetwork$nNodes, rep(.5, self$sampledNetwork$nNodes), directed))
    self$models <- lapply(self$vBlocks, function(i){
      SBM <- switch(paste0(self$family, ifelse(self$sampledNetwork$directed, "Directed", "Undirected")),
                    "BernoulliUndirected" = SBM_BernoulliUndirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                    "BernoulliDirected"   = SBM_BernoulliDirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                    "PoissonUndirected"   = SBM_PoissonUndirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                    "PoissonDirected"     = SBM_PoissonDirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05))
      SBM_VEMfit <- SBM_VEMfit$new(SBM, self$sampledNetwork, self$samplingData)
      SBM_VEMfit$doVEM()
      self$vICLs <- c(self$vICLs, SBM_VEMfit$vICL)
      return(SBM_VEMfit$SBM)
    }
    )
  }
)

SBM_collection$set("public", "smoothing",
                   function() {
                     
                     for(i in rev(seq_along(self$vBlocks))){
                       
                     }
                   }
)

SBM_collection$set("public", "getBestModel",
                   function() {
                     return(nrow(self$models[[which.min(self$vICLs)]]$connectParam))
                   }
)

### Tests :
# SBM :
# mySBM <- SBM_BernoulliUndirected.fit$new(100, c(1/2, 1/2), matrix(c(.5, .05, .05, .5),2,2))
# 
# # # Sampled SBM :
# mySampledSBM   <- sampling_starDegree$new(100, c(-3, .13), FALSE)
# Y              <- mySBM$rSBM()$adjacencyMatrix
# # Z              <- t(rmultinom(100, size = 1, prob = c(.5, .5)))
# # Znum           <- Z %*% c(1:2)
# # sample         <- mySampledSBM$rSampling(Y, Z)
# sample         <- mySampledSBM$rSampling(Y)
# 
# # Sampled SBM 2 (MAR):
# # mySampledSBM  <- sampling_randomPairMAR$new(100, 1/2, FALSE)
# # Y             <- mySBM$rSBM()$adjacencyMatrix
# # sampledNetwork <- mySampledSBM$rSampling(Y)
# 
# # VEM :
# sbm <- SBM_collection$new(sample$adjacencyMatrix, 2, "starDegree", "Bernoulli", FALSE)
# sbmsmoothing()
# 
# sbm$getBestModel()
# sbm$vICLs[-1]









