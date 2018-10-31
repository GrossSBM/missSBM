#' R6 Class definition of an SBM-fit
#'
#' This class is designed to adjust a Stochastic Block Model on a fully observed network.
#' The doVEM method performs inference via Variational EM.
#'
#' This class is virtual: inference is effective only for instance of one of the two child
#' classes 'SBM_fit_nocovariate' and 'SBM_fit_covariates'
#'
#' @import R6
#' @include SBM_fit-Class.R
#' @export
SBM_fit <-
R6::R6Class(classname = "SBM_fit",
  inherit = SBM,
  private = list(
    tau  = NULL  # variational parameters for posterior probablility of class belonging
  ),
  public = list(
    vBound = function(adjMatrix) {self$vExpec(adjMatrix) + self$entropy},
    vICL   = function(adjMatrix) {-2 * self$vExpec(adjMatrix) + self$penalty}
  ),
  active = list(
    blocks      = function(value) {if (missing(value)) return(private$tau) else  private$tau <- value},
    memberships = function(value) {apply(private$tau, 1, which.max)},
    penalty     = function(value) {(self$df_connectParams + self$df_covarParams) * log(self$nDyads) + self$df_mixtureParams * log(self$nNodes)},
    entropy     = function(value) {-sum(xlogx(private$tau))},
    connectProb = function(value) {
      PI <- private$tau %*% private$pi %*% t(private$tau)
      if (self$has_covariates) {
        PI <- logistic(PI + roundProduct(private$phi, private$beta))
      }
      PI
    }
  )
)

SBM_fit$set("public", "doVEM",
  function(adjMatrix, threshold = 1e-4, maxIter = 10, fixPointIter = 3, trace = FALSE) {

    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", maxIter)
    objective <- vector("numeric", maxIter)
    iterate <- 0; stop <- FALSE

    ## Starting the variational EM algorithm
    if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    while (!stop) {
      iterate <- iterate + 1
      if (trace) cat(" iteration #:", iterate, "\r")

      pi_old <- private$pi # save old value of parameters to assess convergence

      # Variational E-Step
      self$update_blocks(adjMatrix, fixPointIter)
      # M-step
      self$update_parameters(adjMatrix)

      # Assess convergence
      delta[iterate] <- sqrt(sum((private$pi - pi_old)^2)) / sqrt(sum((pi_old)^2))
      stop <- (iterate > maxIter) |  (delta[iterate] < threshold)
      objective[iterate] <- self$vBound(adjMatrix)
    }
    if (trace) cat("\n")
    res <- data.frame(delta = delta[1:iterate], objective = objective[1:iterate])
    res
  }
)

