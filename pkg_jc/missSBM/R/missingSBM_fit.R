#' R6 Class definition of an SBM-fit
#'
#' This class is designed to adjust a Stochastic Block Model on a network with missing entries.
#'
#' @include SBM_fit-Class.R
#' @include networkSampling_fit-Class.R
#' @export
missingSBM_fit <-
R6Class(classname = "missingSBM_fit",
  private = list(
    imputedNet = NULL, # adjacency matrix with imputed NA
    sampling   = NULL, # object of class networkSampling_fit
    SBM        = NULL  # object of class SBM_fit
  ),
  public = list(
    initialize = function(adjMatrix = NA, nBlocks = NA, netSampling = NA, clusterInit = "spectral") {

      stopifnot(netSampling %in% available_samplings)

      ## network data with basic imputation
      private$imputedNet <- adjMatrix;
      private$imputedNet[is.na(adjMatrix)] <- mean(private$imputedNet[!is.na(adjMatrix)])

      ## construct the SBM fit
      private$SBM <- SBM_fit$new(private$imputedNet, nBlocks, clusterInit)

      ## construct the network sampling fit
      private$sampling <- switch(netSampling,
        "dyad"            = dyadSampling_fit$new(adjMatrix),
        "node"            = nodeSampling_fit$new(adjMatrix),
        "block"           = blockSampling_fit$new(adjMatrix, private$SBM$blocks),
        "double_standard" = doubleStandardSampling_fit$new(adjMatrix),
        "degree"          = degreeSampling_fit$new(adjMatrix),
        "snowball"        = snowballSampling_fit$new(adjMatrix)
      )
    }
  ),
  active = list(
    imputedNetwork = function(value) {private$imputedNet},
    fittedSBM = function(value) {private$SBM},
    fittedSampling = function(value) {private$sampling}
  )
)
