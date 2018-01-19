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
    sampledNet = NULL, # network data with convenient encoding
    imputedNet = NULL, # imputed network data
    sampling   = NULL, # object of class networkSampling_fit
    SBM        = NULL  # object of class SBM_fit
  ),
  public = list(
    initialize = function(adjMatrix = NA, nBlocks = NA, netSampling = NA, clusterInit = "spectral") {

      stopifnot(netSampling %in% available_samplings)

      ## Turn the adjacency matrix to a proper sampledNetwork object
      private$sampledNet <- sampledNetwork$new(adjMatrix)

      ## network data with basic imputation in NMAR sampling. No imputation need for MAR sampling
      private$imputedNet <- adjMatrix
      # if (netSampling %in% c("double_standard", "degree", "block"))
      #   imputedNet[is.na(adjMatrix)] <- mean(imputedNet[!is.na(adjMatrix)])

      ## construct the SBM fit
      private$SBM <- SBM_fit$new(private$imputedNet, nBlocks, clusterInit)

      ## construct the network sampling fit
      private$sampling <- switch(netSampling,
        "dyad"            = dyadSampling_fit$new(private$sampledNet),
        "node"            = nodeSampling_fit$new(private$sampledNet),
        "block"           = blockSampling_fit$new(private$sampledNet, private$SBM$blocks),
        "double_standard" = doubleStandardSampling_fit$new(private$sampledNet),
        "degree"          = degreeSampling_fit$new(private$sampledNet),
        "snowball"        = snowballSampling_fit$new(private$sampledNet)
      )
    }
  ),
  active = list(
    fittedSBM = function(value) {private$SBM},
    fittedSampling = function(value) {private$sampling}  ,
    sampledNetwork = function(value) {private$sampledNet},
    imputedNetwork = function(value) {private$imputedNet},
    vLogLik = function(value) {private$SBM$vLogLik(private$imputedNet) + private$sampling$vLogLik},
    penalty = function(value) {private$SBM$penalty + private$sampling$penalty},
    vICL = function(value) {-2 * self$vLogLik + self$penalty}
  )
)

#' @export
missingSBM_fit$set("public", "doVEM",
  function(threshold = 1e-5, maxIter = 100, fixPointIter = 5, trace = FALSE) {

    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", maxIter)
    objective <- vector("numeric", maxIter)
    i <- 0; cond <- FALSE
    ## Starting the Variational EM algorithm
    if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    if (trace) cat("\n\tDyads are distributed according to a '",private$SBM$direction,"' SBM with a '" , private$SBM$emissionLaw,"' distribution.\n", sep = "")
    if (trace) cat("\n\tImputation assumes a '", private$sampling$type,"' network-sampling process\n", sep = "")
    while (!cond) {
      i <- i + 1
      if (trace) cat(" iteration #:", i, "\r")

      pi_old <- private$SBM$connectParam # save old value of parameters to assess convergence
      ## ______________________________________________________
      ## Variational E-Step
      #
      # update the variational parameters for block memberships (a.k.a tau)
      private$SBM$update_blocks(private$imputedNet, fixPointIter)
      # update the variational parameters for missing entries (a.k.a nu)
      nu <- private$sampling$update_imputation(private$SBM$blocks, private$SBM$connectParam)
      private$imputedNet[private$sampledNet$NAs] <- nu

      ## ______________________________________________________
      ## M-step
      #
      # update the parameters of the SBM (a.k.a alpha and pi)
      private$SBM$update_parameters(private$imputedNet)
      # update the parameters of network sampling process (a.k.a psi)
      private$sampling$update_parameters(private$imputedNet, private$SBM$blocks)

      ## Check convergence
      delta[i] <- sqrt(sum((private$SBM$connectParam - pi_old)^2)) / sqrt(sum((pi_old)^2))
      cond     <- (i > maxIter) |  (delta[i] < threshold)
      objective[i] <- self$vLogLik
    }
    if (trace) cat("\n")
    res <- data.frame(delta = delta[1:i], objective = objective[1:i])
    res
  }
)
