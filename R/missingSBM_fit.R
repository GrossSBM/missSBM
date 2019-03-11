#' R6 Class definition of an SBM-fit
#'
#' This class is designed to adjust a Stochastic Block Model on a network with missing entries.
#'
#' @include SBM_fit-Class.R
#' @include networkSampling_fit-Class.R
#' @export
missingSBM_fit <-
  R6::R6Class(classname = "missingSBM_fit",
  private = list(
    sampledNet = NULL, # network data with convenient encoding (object of class 'sampledNetwork')
    imputedNet = NULL, # imputed network data (a matrix possibly with NA when MAR sampling is used)
    sampling   = NULL, # fit of the current sampling model (object of class 'networkSampling_fit')
    SBM        = NULL, # fit of the current stochastic block model (object of class 'SBM_fit')
    optStatus  = NULL  # status of the optimization process
  ),
  public = list(
    initialize = function(sampledNet, nBlocks, netSampling,
      clusterInit = ifelse(is.null(covarMatrix), "hierarchical", "spectral"), covarMatrix = NULL, covarSimilarity = l1_similarity) {

      ## Basic arguments checks
      stopifnot(netSampling %in% available_samplings)
      stopifnot(inherits(sampledNet, "sampledNetwork"))
      stopifnot(length(nBlocks) == 1 & nBlocks >= 1 & is.numeric(nBlocks))

      ## Compute the array of covariates, used in all SBM-related computations
      covarArray  <- getCovarArray(covarMatrix, covarSimilarity)

      ## Initial Clustering - Should / Could be a method of sampledNetwork for clarity
      clusterInit <- init_clustering(sampledNet$adjMatrix, nBlocks, covarArray, clusterInit)

      ## network data with basic imputation at startup
      private$imputedNet <- sampledNet$adjMatrix
      Z <- clustering_indicator(clusterInit)
      adjancency0 <- sampledNet$adjMatrix; adjancency0[sampledNet$NAs] <- 0
      pi0 <- check_boundaries((t(Z) %*% adjancency0 %*% Z) / (t(Z) %*% (1 - diag(sampledNet$nNodes)) %*% Z))
      private$imputedNet[sampledNet$NAs] <- (Z %*% pi0 %*% t(Z))[sampledNet$NAs]

      ## Save the sampledNetwork object in the current environment
      private$sampledNet <- sampledNet

      ## Initialize the sampling fit and the SBM fit
      if (is.null(covarMatrix)) {
        private$SBM <- SBM_fit_nocovariate$new(private$imputedNet, clusterInit)
        private$sampling <- switch(netSampling,
          "dyad"            = dyadSampling_fit$new(private$sampledNet),
          "node"            = nodeSampling_fit$new(private$sampledNet),
          "block-node"      = blockSampling_fit$new(private$sampledNet, private$SBM$blocks),
          "double-standard" = doubleStandardSampling_fit$new(private$sampledNet),
          "block-dyad"      = blockDyadSampling_fit$new(private$sampledNet, private$SBM$blocks),
          "degree"          = degreeSampling_fit$new(private$sampledNet, private$SBM$blocks, private$SBM$connectParam)
        )
      } else {
        private$SBM <- SBM_fit_covariates$new(private$imputedNet, clusterInit, covarArray)
        private$sampling <- switch(netSampling,
          "dyad"            = dyadSampling_fit_covariates$new(private$sampledNet, covarArray),
          "node"            = nodeSampling_fit_covariates$new(private$sampledNet, covarMatrix)
        )
      }
    }
  ),
  active = list(
    fittedSBM = function(value) {private$SBM},
    fittedSampling = function(value) {private$sampling}  ,
    sampledNetwork = function(value) {private$sampledNet},
    imputedNetwork = function(value) {private$imputedNet},
    monitoring     = function(value) {private$optStatus},
    entropyImputed = function(value) {
      nu <- private$imputedNet[private$sampledNet$missingDyads]
      res <- -sum(xlogx(nu) + xlogx(1 - nu))
      res
    },
    vBound  = function(value) {private$SBM$vBound(private$imputedNet) + self$entropyImputed + private$sampling$vExpec},
    vExpec  = function(value) {private$SBM$vExpec(private$imputedNet) + private$sampling$vExpec},
    penalty = function(value) {private$SBM$penalty + private$sampling$penalty},
    vICL    = function(value) {-2 * self$vExpec + self$penalty}
  )
)

missingSBM_fit$set("public", "doVEM",
  function(control) {

    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", control$maxIter)
    objective <- vector("numeric", control$maxIter)
    i <- 0; cond <- FALSE
    ## Starting the Variational EM algorithm
    if (control$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    if (control$trace) cat("\n\tDyads are distributed according to a '", private$SBM$direction,"' SBM.\n", sep = "")
    if (control$trace) cat("\n\tImputation assumes a '", private$sampling$type,"' network-sampling process\n", sep = "")
    while (!cond) {
      i <- i + 1
      if (control$trace) cat(" iteration #:", i, "\r")
      pi_old <- private$SBM$connectParam # save current value of the parameters to assess convergence

      ## ______________________________________________________
      ## Variational E-Step
      #
      for (k in seq.int(control$fixPointIter)) {
        # update the variational parameters for missing entries (a.k.a nu)
        nu <- private$sampling$update_imputation(private$SBM$connectProb)
        private$imputedNet[private$sampledNet$NAs] <- nu[private$sampledNet$NAs]
        # update the variational parameters for block memberships (a.k.a tau)
        private$SBM$update_blocks(private$imputedNet, log_lambda = private$sampling$log_lambda)
      }

      ## ______________________________________________________
      ## M-step
      #
      # update the parameters of the SBM (a.k.a alpha and pi)
      private$SBM$update_parameters(private$imputedNet)
      # update the parameters of network sampling process (a.k.a psi)
      private$sampling$update_parameters(private$imputedNet, private$SBM$blocks)

      ## Check convergence
      delta[i] <- sqrt(sum((private$SBM$connectParam - pi_old)^2)) / sqrt(sum((pi_old)^2))
      cond     <- (i > control$maxIter) |  (delta[i] < control$threshold)
      objective[i] <- self$vBound

    }
    if (control$trace) cat("\n")
    private$optStatus <- data.frame(delta = delta[1:i], objective = c(NA, objective[2:i]), iteration = 1:i)
    invisible(private$optStatus)
  }
)

#' @export
missingSBM_fit$set("public", "show",
function() {
  cat("missSBM-fit\n")
  cat("==================================================================\n")
  cat("Structure for storing a SBM fitted under missing data condition   \n")
  cat("==================================================================\n")
  cat("* Useful fields (most are special obejct themselves)              \n")
  cat("  $fittedSBM (the adjusted stochastoc bloc model)                 \n")
  cat("  $fittedSampling (the estimated sampling process)                \n")
  cat("  $sampledNetwork (the sampled network data)                      \n")
  cat("  $imputedNetwork (the adjacency matrix with imputed values)      \n")
  cat("  $monitoring, $vICL, $vBound, $vExpec, $penalty                  \n")
})
missingSBM_fit$set("public", "print", function() self$show())
