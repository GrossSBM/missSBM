#' A collection of adjusted Stochastic Block Model
#'
#' @export
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
    require(parallel)
    self$family   <- family
    self$vBlocks  <- vBlocks
    self$sampling <- sampling
    self$sampledNetwork <- sampledNetwork$new(sample, directed)
    self$samplingData   <- switch(sampling,
                            "doubleStandard" = sampling_doubleStandard$new(self$sampledNetwork$nNodes, c(.5,.5), directed),
                            "class"          = sampling_class$new(self$sampledNetwork$nNodes, .5, directed),
                            "starDegree"     = sampling_starDegree$new(self$sampledNetwork$nNodes, coefficients(glm(self$sampledNetwork$samplingVector~rowSums(sample, na.rm=TRUE), family = binomial(link = "logit"))), directed),
                            "MAREdge"        = sampling_randomPairMAR$new(self$sampledNetwork$nNodes, .5, directed),
                            "MARNode"        = sampling_randomNodesMAR$new(self$sampledNetwork$nNodes, .5, directed),
                            "snowball"       = sampling_snowball$new(self$sampledNetwork$nNodes, rep(.5, self$sampledNetwork$nNodes), directed))
    self$models <- mclapply(self$vBlocks, function(i){
      SBM <- switch(paste0(self$family, ifelse(self$sampledNetwork$directed, "Directed", "Undirected")),
                    "BernoulliUndirected" = SBM_BernoulliUndirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                    "BernoulliDirected"   = SBM_BernoulliDirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                    "PoissonUndirected"   = SBM_PoissonUndirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05),
                    "PoissonDirected"     = SBM_PoissonDirected.fit$new(self$sampledNetwork$nNodes, rep(1, i)/i, diag(.45,i)+.05))
      if(self$family == "Poisson") { init <- "SpectralC" }
      SBM_VEMfit <- SBM_VEMfit$new(SBM, self$sampledNetwork, self$samplingData, init)
      if(self$family == "Bernoulli"){
        SBM_VEMfit$doVEM()
      } else {
        SBM_VEMfit$doVEMPoisson()
      }
      return(SBM_VEMfit)
    }
    , mc.cores = 5)
    self$vICLs <- sapply(self$models, function(x){x$vICL})
  }
)

SBM_collection$set("public", "smoothingBackward",
                   function() {
                     for(i in rev(self$vBlocks[-1])){
                       comb <- combn(i, 2, simplify = FALSE)
                       for(j in 1:length(comb)){
                         cl_fusion <- factor(apply(self$models[[i]]$blockVarParam, 1, which.max))
                         if(length(levels(cl_fusion)) == i){
                           levels(cl_fusion)[which(levels(cl_fusion) == paste(comb[[j]][1]))] <- paste(comb[[j]][2])
                           levels(cl_fusion) <- as.character(1:(i-1))
                           clone             <- self$models[[i-1]]$clone()
                           clone$blockInit   <- as.numeric(cl_fusion)
                           clone$doVEM()

                           if(clone$vICL < self$models[[i-1]]$vICL){
                             self$models[[i-1]] <- clone
                             cat('+')
                           }
                         }
                       }
                     }
                     self$vICLs <- sapply(self$models, function(x){x$vICL})
                   }
)

# SBM_collection$set("public", "smoothingForward",
#                    function() {
#                      for(i in self$vBlocks[-length(self$vBlocks)]]){
#                        for(j in 1:i){
#                        }
#                      }
#                      self$vICLs <- sapply(self$models, function(x){x$vICL})
#                    }
# )

SBM_collection$set("private", "CAH",
                   function(X, Q) {
                     require(hclust)
                     D  <- as.matrix(dist(X, method="manhattan"))
                     D[which(X == 1)] <- D[which(X == 1)] - 2
                     is (anyNA(D))
                       D <- as.dist(additive(D))
                     return(cutree(hclust(as.dist(D), method="ward.D"), Q))
                   }
)



SBM_collection$set("public", "getBestModel",
                   function() {
                     return(nrow(self$models[[which.min(self$vICLs)]]$SBM$connectParam))
                   }
)

