#' an SBM model
#'
#' @field 
#' @field 
#' @field 
#' @field 
#'
#' @importFrom R6 R6Class
#' @export
SBM_BernoulliUndirectedfit <-
  R6Class(classname = "SBM_BernoulliUndirectedfit",
          inherit = SBM_BernoulliUndirected,
          public = list(
            ## fields

            ## methods
            initialize = function(nNodes=NA, mixtureParam=NA, connectParam=NA) {
              self$nNodes       <- nNodes
              self$mixtureParam <- mixtureParam
              self$connectParam <- connectParam
              self$nBlocks      <- length(mixtureParam)
            }
          )
  )