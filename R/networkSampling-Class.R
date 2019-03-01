#' definition of R6 Class 'networkSampling'
#'
#' this virtual class is the mother of all subtypes of networkSampling (either sampler or fit)
#'
#' @include utils.R
#'
#' @import R6
networkSampling <-
R6::R6Class(classname = "networkSampling",
  ## fields
  private = list(
    name  = NULL, # type of sampling
    psi   = NULL  # vector of missing parameters
  ),
  public = list(
    ## methods
    initialize = function(type = NA) {
      stopifnot(type %in% available_samplings)
      private$name <- type
    },
    show = function(model = paste0(private$name, "-model for network sampling\n")) {
      cat(model)
      cat("==================================================================\n")
      cat("Structure for handling network sampling in missSBM.\n")
      cat("==================================================================\n")
      cat("* Useful fields \n")
      cat("  $type, $parameters, $df\n")
    },
    print = function() {self$show()}
  ),
  active = list(
    type = function(value) {private$name},
    parameters = function(value) {private$psi},
    ## degree of freedom are just the size of the vector of missing parameters
    df = function(value) {if(is.matrix(private$psi)){return(length(private$psi[lower.tri(private$psi, diag = TRUE)]))} else {length(private$psi)}}
  )
)

#' @include utils.R
#'
#' @import R6
networkSamplingCovariates <-
R6::R6Class(classname = "networkSamplingCovariates",
  inherit = networkSampling,
  ## fields
  private = list(
    X   = NULL,
    phi = NULL # array of covariates
  ),
  public = list(
    ## methods
    initialize = function(type = NA, covariates = NA) {
      stopifnot(type %in% available_samplings_covariates)
      private$X <- covariates
      N <- nrow(covariates)
      M <- ncol(covariates)
      phi <- array(dim = c(N, N, M))
      for (i in 1:N)
        for (j in 1:N)
          phi[i,j,] <- -abs(covariates[i, ] - covariates[j, ])
      private$phi <- phi
      private$name <- type
    }
  )
)
