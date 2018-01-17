#' definition of R6 Class 'sampling'
#'
#' this virtual class is the mother of all subtypes of sampling (either model or fit)
#'
#' @include utils.R
#'
#' @import R6
#' @export
sampling <-
R6Class(classname = "sampling",
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
    }
  ),
  active = list(
    type = function(value) {
      if (missing(value)) return(private$psi) else private$psi <- value
    },
    parameters = function(value) {
      if (missing(value)) return(private$psi) else private$psi <- value
    },
    ## degree of freedom are just the size of the vector of missing parameters
    df = function(value) {length(private$psi)}
  )
)
