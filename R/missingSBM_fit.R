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
    sampledNet = NULL, # network data with convenient encoding (object of class 'sampledNetwork')
    imputedNet = NULL, # imputed network data (a matrix possibly with NA when MAR sampling is used)
    sampling   = NULL, # fit of the current sampling model (object of class 'networkSampling_fit')
    SBM        = NULL  # fit of the current stochastic block model (object of class 'SBM_fit')
  ),
  public = list(
    initialize = function(sampledNet, nBlocks, netSampling, clusterInit = "spectral") {

      ## Basic arguments checks
      stopifnot(netSampling %in% available_samplings)
      stopifnot(inherits(sampledNet, "sampledNetwork"))
      stopifnot(length(nBlocks) == 1 & nBlocks >= 1 & is.numeric(nBlocks))

      ## Save the sampledNetwork object in the current environment
      private$sampledNet <- sampledNet

      ## network data without any imputation at startup
      private$imputedNet <- sampledNet$adjacencyMatrix

      ## construct and initialize the SBM fit
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
    entropyImputed = function(value) {
      nu <- private$imputedNet[private$sampledNet$NAs]
      res <- -sum(xlogx(nu) + xlogx(1 - nu))
      res
    },
    vBound  = function(value) {private$SBM$vBound(private$imputedNet) + self$entropyImputed + private$sampling$logLik},
    penalty = function(value) {private$SBM$penalty + private$sampling$penalty},
    vBIC = function(value) {-2 * self$vBound + self$penalty},
    vICL = function(value) {-2 * (private$SBM$vExpec(private$imputedNet) + private$sampling$logLik) + self$penalty}
  )
)

#' @export
missingSBM_fit$set("public", "doVEM",
  function(threshold = 1e-4, maxIter = 200, fixPointIter = 5, trace = FALSE) {

    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", maxIter)
    objective <- vector("numeric", maxIter)
    i <- 0; cond <- FALSE
    ## Starting the Variational EM algorithm
    if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    if (trace) cat("\n\tDyads are distributed according to a '", private$SBM$direction,"' SBM with a '" , private$SBM$emissionLaw,"' distribution.\n", sep = "")
    if (trace) cat("\n\tImputation assumes a '", private$sampling$type,"' network-sampling process\n", sep = "")
    while (!cond) {
      i <- i + 1
      if (trace) cat(" iteration #:", i, "\r")

      pi_old <- private$SBM$connectParam # save current value of the parameters to assess convergence
      ## ______________________________________________________
      ## Variational E-Step
      #
      # update the variational parameters for missing entries (a.k.a nu)
      nu <- private$sampling$update_imputation(private$SBM$blocks, private$SBM$connectParam)
      private$imputedNet[private$sampledNet$NAs] <- nu[private$sampledNet$NAs]
      # update the variational parameters for block memberships (a.k.a tau)
      private$SBM$update_blocks(private$imputedNet, fixPointIter, log_lambda = private$sampling$log_lambda)

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
      objective[i] <- self$vBound
    }
    if (trace) cat("\n")
    res <- data.frame(delta = delta[1:i], objective = c(NA, objective[2:i]))
    res
  }
)

#' @import igraph
#' @export
missingSBM_fit$set("public", "plot",
  function(type = c("imputedNetwork", "connectivity")) {
    type <- match.arg(type)
    if (type == "imputedNetwork") {
      gg_image_NA(private$imputedNet, private$SBM$memberships)
    }
    if (type == "connectivity") {
      plot(
        graph_from_adjacency_matrix(
          private$SBM$connectParam,
          mode = ifelse(private$sampledNet$is_directed, "directed", "undirected"),
          weighted = TRUE, diag = TRUE
        ), main = "SBM connectivity matrix"
      )
    }
  }
)
