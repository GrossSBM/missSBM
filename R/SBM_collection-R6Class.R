#' a collection of adjusted Stochastic Block Model
#'
#' @importFrom R6 R6Class
SBM_collection <-
  R6Class(classname = "SBM_collection",
    public = list(
      sampledNetwork = NULL, # the sampled network data
      vBlocks  = NULL, # the vector of number of blocks considered in the collection
      models   = NULL, # the collection of SBM fit
      sampling = NULL # the sampling design for missing data modeling
      ## methods
      ### TODO: VE-step and M-step
      ## criteria = NULL, # variational bound and ICL for all models
  )
)
#
# SBM_collection$set("public", "doVEstep",
#   function(model) {
#
#   }
# )
#
# SBM_collection$set("public", "doMstep",
#   function(model) {
#
#   }
# )
#
