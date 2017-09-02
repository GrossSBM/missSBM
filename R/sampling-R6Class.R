#' an SBM model
#'
#' @field nNodes          number of nodes
#' @field nBlocks         number of blocks
#' @field blockProportion vector of block proportion (a.k.a. alpha)
#' @field modelParameters vector of model parameters (a.k.a. theta)
#'
#' @importFrom R6 R6Class
#' @export
sampling <-
R6Class(classname = "sampling",
  public = list(
    ## fields
    nNodes         = NULL, # number of nodes
    missingParam   = NULL, # vector of missing parameters (a.k.a. alpha)
    completeLogLik = NULL, #
    ## methods
    initialize = function(nNodes=NA, missingParam=NA) {
      self$nNodes       <- nNodes
      self$missingParam <- missingParam
    }
  )
)

### TODO
#
# sampling_doubleStandard <-
# R6Class(classname = "sampling_doublestandard",
#   public = list(
#    initialize = function() {},
#    rSampling = function(_some parameters_) {},
#    samplingLogLik = function() {},
#   )
# )
#
# TODO : idem pour starDegree, blockCentered
