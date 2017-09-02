#' an SBM fit, i.e. an adjusted SBM
#'
#' @field completedNetwork
#' @field missingDyadsProb
#' @field blocksProb
#' @field lowerBound
#' @field ICL
#'
#' @include SBM-R6class.R
#'
#' @importFrom R6 R6Class
#' @export
SBM_VEMfit <-
R6Class(classname = "SBM_VEMfit",
  public = list(
    ## fields
    completedNetwork = NULL, # the completed adjacency matrix of the initial network
    missingVarParam  = NULL, # variational parameters for missing entries (a.k.a. nu)
    blockVarParam    = NULL, # variational parameters for latent blocks, (a.k.a. tau)
    SBM              = NULL, #
    ## methods
    lowerBound       = NULL, # variational lower bound (a.k.a. J)
    vICL             = NULL, # compute the (variational) integrated complete likelihood
    blocks           = NULL  # get the most probable blocks
  )
)

