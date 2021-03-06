#' Definition of R6 Class 'networkSampling'
#'
#' this virtual class is the mother of all subtypes of networkSampling (either sampler or fit)
#' It is used to define a sampling model for a network.
#' It has a rSampling method which takes an adjacency matrix as an input and send back an object with class partlyObservedNetwork.
#'
#' @include utils_missSBM.R
#'
#' @import R6
networkSampling <-
R6::R6Class(classname = "networkSampling",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    name  = NULL, # type of sampling process
    psi   = NULL, # vector of parameters of the process originating the missing values
    rho   = NULL  # the probability for sampling either a dyad or a node
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling
    #' @param type character for the type of sampling. must be in ("dyad", "covar-dyad", "node", "covar-node", "block-node", "block-dyad", "double-standard", "degree")
    #' @param parameters the vector of parameters associated to the sampling at play
    initialize = function(type = NA, parameters = NA) {
      stopifnot(type %in% available_samplings)
      private$name <- type
      private$psi  <- parameters
    },
    #' @description show method
    #' @param type character used to specify the type of sampling
    show = function(type = paste0(private$name, "-model for network sampling\n")) {
      cat(type)
      cat("==================================================================\n")
      cat("Structure for handling network sampling in missSBM.\n")
      cat("==================================================================\n")
      cat("* Useful fields \n")
      cat("  $type, $parameters, $df\n")
    },
    #' @description User friendly print method
    print = function() { self$show() }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field type a character for the type of sampling
    type = function(value) {private$name},
    #' @field parameters the vector of parameters associated with the sampling at play
    parameters = function(value) {private$psi},
    #' @field df the number of entries in the vector of parameters
    df = function(value) {length(private$psi)}
  )
)
