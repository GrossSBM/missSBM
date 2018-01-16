#' @include sampling_model-Class.R
#' @import R6
#' @export
sampling_fit <-
R6Class(classname = "sampling_fit",
  inherit = sampling_model,
  public = list(
    initialize = function(adjMatrix) {
      private$net <- sampledNetwork$new(adjMatrix)
    }
  ),
  active = list(
    ## degree of freedom are just the size of the vector of missing parameters
    df = function(value) {length(private$psi)},
    ## nDyads automatically handle the directed/undirected cases
    penalty = function(value) {log(private$net$nDyads) * self$df}
  )
)

#' @export
sampling_fit_dyad <-
R6Class(classname = "sampling_fit_dyad",
  inherit = sampling_fit,
  private = list(
    card_D_o = NULL, # stats required by the likelihood
    card_D_m = NULL  # number of observed, repectively missing dyads
  ),
  public = list(
    initialize = function(adjMatrix) {
      super$initialize(adjMatrix)
      private$name <- "dyad"
      private$card_D_o <- length(private$net$observedDyads)
      private$card_D_m <- length(private$net$missingDyads )
      private$psi <- private$card_D_o / private$net$nDyads
    }
  ),
  active = list(
    vLogLik = function() {
      res <- (private$card_D_o - private$net$nNodes) * logx(private$psi) + private$card_D_m * log1mx(private$psi)
      res
    }
  )
)

#' @export
sampling_fit_node <-
R6Class(classname = "sampling_fit_node",
  inherit = sampling_fit,
  private = list(
    card_N_o = NULL, # stats required by the likelihood
    card_N_m = NULL  # number of observed, repectively missing nodes
  ),
  public = list(
    initialize = function(adjMatrix) {
      super$initialize(adjMatrix)
      private$name <- "node"
      private$card_N_o <- sum( private$net$observedNodes)
      private$card_N_m <- sum(!private$net$observedNodes)
      private$psi <- private$card_N_o / (private$card_N_o + private$card_N_m)
    }
  ),
  active = list(
    vLogLik = function() {
### SHOULD'NT BE A SQUARE FOR N_obs SOMEWHERE?
### THE LOGLIK IS HIGH COMPARED TO THE ONE IN THE DYAD CASE...
      res <- private$card_N_o * logx(private$psi) + private$card_N_m * log1mx(private$psi)
      res
    }
  )
)

#' @export
sampling_fit_double_standard <-
R6Class(classname = "sampling_fit_double_standard",
  inherit = sampling_fit,
  private = list(
    So     = NULL, ## statistics only requiring the observed part of the network
    So.bar = NULL, ## can be computed once for all during the initialization
    Sm     = NULL, ## these ones will be updated during the algorithm
    Sm.bar = NULL
  ),
  public = list(
    initialize = function(adjMatrix, missingInit = NA) {
      super$initialize(adjMatrix)
      private$name <- "double_standard"
      private$So     <- sum(    private$net$adjacencyMatrix[private$net$observedDyads])
      private$So.bar <- sum(1 - private$net$adjacencyMatrix[private$net$observedDyads])
      ## SEE HOW TO "COMPLETE" THE NETWORK AT START-UP IN ORDER TO INITIALIZE PSI
      if (is.na(missingInit))
        missingInit <- rep(mean(private$net$adjacencyMatrix, na.rm = TRUE), length(private$net$missingDyads))
      self$update_missing(missingInit)
    },
    update_missing = function(nu) {
      private$Sm     <- sum(    nu)
      private$Sm.bar <- sum(1 - nu)
      private$psi    <- c(private$So.bar / (private$So.bar + private$Sm.bar), private$So / (private$So + private$Sm))
    }
  ),
  active = list(
    vLogLik = function(value) {
      res <- logx(private$psi[2]) * private$So + logx(private$psi[1]) * private$So.bar +
        log1mx(private$psi[2]) * private$Sm + log1mx(private$psi[1]) * private$Sm.bar
      res
    }
  )
)

#' @export
sampling_fit_block <-
R6Class(classname = "sampling_fit_block",
  inherit = sampling_fit,
  private = list(
    So = NULL, ## sum_(i in Nobs ) Z_iq
    Sm = NULL  ## sum_(i in Nmiss) Z_iq
  ),
  public = list(
    initialize = function(adjMatrix, blockInit) {
      super$initialize(adjMatrix)
      private$name <- "block"
      self$update_missing(blockInit)
    },
    update_missing = function(Z) {
      private$So <- colSums(Z[ private$net$observedNodes, ])
      private$Sm <- colSums(Z[!private$net$observedNodes, ])
      private$psi <- private$So / (private$So + private$Sm)
    }
  ),
  active = list(
    vLogLik = function() {
      res <- c(crossprod(private$So, log(private$psi)) +  crossprod(private$Sm, log(1 - private$psi)))
      res
    }
  )
)

sampling_fit_degree <-
R6Class(classname = "sampling_fit_degree",
  inherit = sampling,
  public = list(
    vLogLik = function(completedNetwork) {
      prob <- logistic(private$psi[1] + private$psi[2] * rowSums(completedNetwork))
      res  <- log( prob^private$net$observedNodes %*% (1 - prob)^(!private$net$observedNodes) )
      res
    }
    # updatePsi = function(completedNetwork, sampledNetwork, blockVarParam, taylorVarParam) {
    #   networkWithZeros     <- completedNetwork
    #   networkWithZeros[sampledNetwork$missingDyads] <- 0
    #   Dtilde <- rowSums(completedNetwork)
    #   Dchap  <- rowSums((completedNetwork-networkWithZeros)*(1-(completedNetwork-networkWithZeros))) + Dtilde^2
    #   Nmiss  <- length(which(sampledNetwork$samplingVector == 0))
    #   b1     <- ( (((2*sum(private$g(taylorVarParam)*Dtilde))*(-length(Nmiss) + 0.5*sampledNetwork$nNodes)))/(sum(private$g(taylorVarParam))) - (-sum(Dtilde[Nmiss]) + sum(Dtilde)*0.5) )
    #   b2     <- ( 2*sum(private$g(taylorVarParam)*Dchap) - (((2*sum(private$g(taylorVarParam)*Dtilde))^2 ))/(sum(private$g(taylorVarParam))))
    #   b      <- b1/b2
    #   a      <- -(b*(2*sum(private$g(taylorVarParam)*Dtilde)) + (-length(Nmiss) + 0.5*sampledNetwork$nNodes))/(sum(private$g(taylorVarParam)))
    #   psi    <- c(a,b)
    #   return(psi)
    # },
    # penality = function(nBlocks) {
    #   if(self$directed){
    #     return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + 2*(nBlocks-1)*log(self$nNodes))
    #   } else {
    #     return((nBlocks*(nBlocks+1)/2)*log(self$nNodes*(self$nNodes-1)/2) + 2*(nBlocks-1)*log(self$nNodes))
    #   }
    # }
  )
)

sampling_fit_snowball <-
R6Class(classname = "sampling_fit_snowball",
  inherit = sampling,
  public = list(
    samplingLogLik = function() {
      0
    },
    # samplingLogLik = function(sampledNetwork, completedNetwork) {
    #   return(log((self$missingParam^sampledNetwork$samplingVector)%*%((1-self$missingParam)^(1-sampledNetwork$samplingVector))))
    # },
    penality = function(nBlocks) {
      if(self$directed){
        return((nBlocks^2)*log(self$nNodes*(self$nNodes-1)) + nBlocks*log(self$nNodes))
      }
      else {
        return(nBlocks*(nBlocks+1)/2*log(self$nNodes*(self$nNodes-1)/2) + nBlocks*log(self$nNodes))
      }
    }
  )
)

