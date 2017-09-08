#' an SBM fit, i.e. an adjusted SBM
#'
#' @field completedNetwork
#' @field missingDyadsProb
#' @field blocksProb
#' @field lowerBound
#' @field ICL
#'
#' @include SBM-R6class.R
#'
#' @importFrom R6 R6Class
#' @export
SBM_VEMfit <-
  R6Class(classname = "SBM_VEMfit",
          public = list(
            ## fields
            completedNetwork = NULL, # the completed adjacency matrix of the initial network
            missingVarParam  = NULL, # variational parameters for missing entries (a.k.a. nu)
            blockVarParam    = NULL, # variational parameters for latent blocks, (a.k.a. tau)
            SBM              = NULL, #
            sampling         = NULL, #
            controlVEstep    = NULL,
            controlMstep     = NULL,
            ## methods
            lowerBound       = NULL, # variational lower bound (a.k.a. J)
            vICL             = NULL, # compute the (variational) integrated complete likelihood
            blocks           = NULL,  # get the most probable blocks
            fixPoint         = NULL,
            updateNu         = NULL,
            VEstep           = function(completedNetwork) {
              self$blockVarParam <- self$fixPoint(completedNetwork, SBM, self$blockVarParam)
              if(!(class(sampling) %in% c("sampling_randomPairMAR", "sampling_randomNodesMAR", "sampling_snowball"))){
                self$missingVarParam <- updateNu(SBM, sampling, self$blockVarParam)
              }
            },
            Mstep            = 
            
            
            initialize = function(SBM, Blocks) {

              #### à compléter

              self$SBM <- SBM
              self$blockVarParam <- matrix(0,SBM$nNodes,SBM$nBlocks)
              self$blockVarParam[cbind(1:SBM$nNodes, Blocks)] <- 1

              #### à compléter
            }
          )
  )



# doVE_step = function(completedNetwork) {
#
#   ## E-step: fixed point to solve for Tau
#   cond.FP <- FALSE; iter <- 0
#   while(!cond.FP) {
#     iter <- iter + 1
#     if(is.null(CL)){
#
#       completedNetwork %*% self$blockVarParam %*% t(log(self$SBM$connecParam)) + X3 %*% self$blockVarParam %*% t(log(1-self$SBM$connecParam))
#
#       self$blockVarParam <- exp(sweep(X2 %*% self$blockVarParam %*% t(log(pi)) + X3 %*% self$blockVarParam %*% t(log(1-pi)),2,log(alpha),"+"))
#       
#       num <- rowSums(Tau)
#       Tau <- Tau/num
#     }
#     Tau[is.nan(Tau)] <- 0.5
#     cond.FP <- (iter > maxIter.FP) | (sum((Tau.old - Tau)^2)/sum(Tau^2) < eps.FP)
#
#     ## MaJ des nu_ij
#     PI       <- log(pi) - log(1-pi)
#     X2[miss] <- (logistic( log(1-psi[2]) - log(1-psi[1]) + Tau %*% PI %*% t(Tau)))[miss]
#     X3       <- 1-X2; diag(X3) <- 0
#   }
#
#   ## Storing the current estimates
#   theta[[i]]   <- list(pi=pi, alpha=alpha)
#   psi.all[[i]] <- psi
#   Tau.all[[i]] <- Tau
#
#
# },