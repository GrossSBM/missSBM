#' Definition of R6 Class 'networkSampling'
#'
#' this virtual class is the mother of all subtypes of networkSampling (either sampler or fit)
#' It is used to define a sampling model for a network.
#' It has a rSampling method which takes an adjacency matrix as an input and send back an object with class sampledNetwork.
#'
#' @include utils_missSBM.R
#'
#' @import R6
networkSampling <-
R6::R6Class(classname = "networkSampling",
  ## fields
  private = list(
    name  = NULL, # type of sampling
    psi   = NULL, # vector of missing parameters
    rho   = NULL  # the prior probability for sampling either and dyad or an node
  ),
  public = list(
    ## methods
    initialize = function(type = NA, parameters = NA) {
      stopifnot(type %in% available_samplings)
      private$name <- type
      private$psi  <- parameters
    },
    show = function(model = paste0(private$name, "-model for network sampling\n")) {
      cat(model)
      cat("==================================================================\n")
      cat("Structure for handling network sampling in missSBM.\n")
      cat("==================================================================\n")
      cat("* Useful fields \n")
      cat("  $type, $parameters, $prob, $df\n")
    },
    print = function() {self$show()}
  ),
  active = list(
    type = function(value) {private$name},
    prob  = function(value) {private$rho},
    parameters = function(value) {private$psi},
    df = function(value) {length(private$psi)}
  )
)
