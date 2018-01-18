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
    sampling   = NULL, # object of class networkSampling_fit
    SBM        = NULL  # object of class SBM_fit
  ),
  public = list(
    initialize = function(adjMatrix = NA, nBlocks = NA, netSampling = NA, clusterInit = "spectral") {

      stopifnot(netSampling %in% available_samplings)

      ## network data with basic imputation
      imputedNet <- adjMatrix;
      imputedNet[is.na(adjMatrix)] <- mean(imputedNet[!is.na(adjMatrix)])

      ## construct the SBM fit
      private$SBM <- SBM_fit$new(imputedNet, nBlocks, clusterInit)

      ## construct the network sampling fit
      private$sampling <- switch(netSampling,
        "dyad"            = dyadSampling_fit$new(adjMatrix),
        "node"            = nodeSampling_fit$new(adjMatrix),
        "block"           = blockSampling_fit$new(adjMatrix, private$SBM$blocks),
        "double_standard" = doubleStandardSampling_fit$new(adjMatrix),
        "degree"          = degreeSampling_fit$new(adjMatrix),
        "snowball"        = snowballSampling_fit$new(adjMatrix)
      )
    },
    vLogLik = function() {
      private$fittedSBM$vLogLik(imputedNet)
    }
  ),
  active = list(
    fittedSBM = function(value) {private$SBM},
    fittedSampling = function(value) {private$sampling}
  )
)

#' @export
missingSBM_fit$set("public", "doVEM",
  function(threshold = 1e-5, maxIter = 100, fixPointIter = 5, trace = FALSE) {

    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", maxIter); delta[1] <- NA
    objective <- vector("numeric", maxIter)
    i <- 0; cond <- FALSE
    ## Starting the Variational EM algorithm
    if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    if (trace) cat("\n Missing entries are imputed based on the assumption of a ",private$fittedSampling$type," network sampling process\n")
    while (!cond) {
      i <- i + 1
      if (trace) cat("iteration #:", i, "\r")

      pi_old <- private$SBM$connectParam # save old value of parameters to assess convergence

      ## ______________________________________________________
      ## M-step

      # update the parameters of the SBM (a.k.a alpha and pi)
      private$SBM$update_parameters(private$sampling$imputedNetwork)
      # update the parameters of network sampling process (a.k.a psi)
      private$sampling$update_parameters(private$SBM$blocks)

      ## ______________________________________________________
      ## Variational E-Step

      # update the variational parameters for block memberships (a.k.a tau)
      private$SBM$update_blocks(private$sampling$imputedNetwork, fixPointIter)
      # update the variational parameters for missing entries (a.k.a nu)
      private$sampling$update_imputation(private$SBM$blocks, private$SBM$connectParam)

      ## Check convergence
      if (i > 1) {
        delta[i] <- sqrt(sum((private$SBM$connectParam - pi_old)^2)) / sqrt(sum((pi_old)^2))
        cond     <- (i > maxIter) |  (delta[i] < threshold)
      }
      objective[i] <- NA
    }
    if (trace) cat("\n")
    res <- data.frame(delta = delta[1:i], objective = objective[1:i])
    res
  }
)
