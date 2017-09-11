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
            sampledNetwork   = NULL, # 
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
            VEstep           = function() {
              self$blockVarParam <- self$fixPoint(self$SBM, self$blockVarParam, self$completedNetwork)
              if(!(class(self$sampling) %in% c("sampling_randomPairMAR", "sampling_randomNodesMAR", "sampling_snowball"))){
                self$missingVarParam <- updateNu(self$SBM, self$sampling, self$blockVarParam)
              }
            },
            Mstep            = function() {
              if(!(class(self$sampling) %in% c("sampling_randomPairMAR", "sampling_randomNodesMAR", "sampling_snowball"))){
                self$SBM$connectParam <- self$maximization(self$SBM, self$completedNetwork, self$blockVarParam)
              } else {
                self$maximization(self$SBM, self$completedNetwork, self$blockVarParam, self$sampledNetwork$samplingMatrix)
              }
            },
            initialize = function(SBM, Blocks) {
              
              
              self$SBM <- SBM
              self$blockVarParam <- matrix(0,SBM$nNodes,SBM$nBlocks)
              self$blockVarParam[cbind(1:SBM$nNodes, Blocks)] <- 1
              
            }
          )
  )

SBM_VEMfit$set("public", "estimate",
               doVEM = function(completedNetwork) {
                 
                 n <- nrow(X)
                 miss <- which(is.na(X))
                 obs  <- which(!is.na(X))
                 ## A. INITIALISATION DES PARAMÈTRES
                 
                 ## A.1 PARAMÈTRES CONTRÔLANT L'OPTIMISATION
                 eps        <- 1e-5
                 maxIter    <- 250
                 eps.FP     <- 1e-2
                 maxIter.FP <- 5
                 mc.cores   <- mc.cores
                 zero       <- .Machine$double.eps
                 if(is.character(cl.init)) stopifnot(cl.init %in% c("mar.CAH","CAH", "spectral", "spectral.CAH"))
                 
                 ## A.2 PARAMÈTRES FIXES EN COURS D'ALGORITHME (données)
                 R  <- 1*(!is.na(X)); diag(R) <- 1         # matrix of observed data
                 X1 <- X; X1[is.na(X)] <- 0                # indicator for truly observed edges
                 # X0 <- (1-X1); diag(X0) <- 0; X0 <- X0 * R # indicator for truly unobserved edges
                 X2 <- X
                 
                 models <- mclapply(seq.Q, function(Q) {
                   ## ===================================================================
                   ## B. VEM À NOMBRE DE CLASSES Q FIXÉ
                   
                   ## B.1 INITIALISATION DES CLUSTERS
                   if(is.null(CL)){
                     if (is.character(cl.init) | is.null(cl.init)) { ## CAH, spectral clustering or user defined
                       cl0 <- switch(cl.init,
                                     "spectral"     = SpectralClustering(X, Q),
                                     "CAH"          = graphCAH(X, Q),
                                     "mar.spectral" = SpectralClustering(X, Q),
                                     "mar.CAH"      = graphCAH(X, Q))
                     } else {
                       if (is.list(cl.init)) {
                         cl0 <- cl.init[[Q]]
                       } else {
                         cl0 <- cl.init[[which(seq.Q == Q)]]
                       }
                     }
                   } else {
                     cl0 <- CL
                   }
                   if(!is.null(CL.init) & is.null(CL)){
                     cl0 <- CL.init
                   }
                   
                   Tau      <- matrix(0,n,Q); Tau[cbind(1:n, cl0)] <- 1
                   piInitNu <- (t(Tau)%*% X1 %*%Tau) / (t(Tau)%*%((1-diag(n)))%*%Tau)
                   piInitNu[is.nan(piInitNu)] <- zero ; piInitNu[piInitNu > 1-zero] <- 1-zero ; piInitNu[piInitNu < zero] <- zero
                   X2[miss] <- ((Tau) %*% piInitNu %*% t(Tau))[miss]
                   X3       <- 1-X2; diag(X3) <- 0
                   
                   theta   <- vector("list", length = maxIter)
                   psi.all <- vector("list", length = maxIter)
                   Tau.all <- vector("list", length = maxIter)
                   J       <- vector("numeric", maxIter)
                   conv    <- vector("numeric", maxIter); conv[1] <- NA
                   
                   
                   ## B.2 EXPECTATION-MAXIMIZATION ALGORITHM
                   i <- 0; cond <- FALSE
                   while(!cond){
                     i <- i+1
                     ## M-step explicit form of pi and alpha
                     pi             <- (t(Tau)%*% X2 %*%Tau) / (t(Tau)%*%((1-diag(n)))%*%Tau)
                     pi[is.nan(pi)] <- zero ; pi[pi > 1-zero] <- 1-zero ; pi[pi < zero] <- zero
                     alpha          <-  colMeans(Tau)
                     
                     SS  <- c(sum(1-X2[obs])-n,sum(X2[obs]))
                     aux <- c(sum(1-X2     )-n,sum(X2)     )
                     psi <- SS/aux
                     
                     ## E-step: fixed point to solve for Tau
                     cond.FP <- FALSE; iter <- 0
                     while(!cond.FP){
                       Tau.old <- Tau
                       iter <- iter + 1
                       if(is.null(CL)){
                         Tau <- exp(sweep(X2 %*% Tau %*% t(log(pi)) + X3 %*% Tau %*% t(log(1-pi)),2,log(alpha),"+"))
                         num <- rowSums(Tau)
                         Tau <- Tau/num
                       }
                       Tau[is.nan(Tau)] <- 0.5
                       cond.FP <- (iter > maxIter.FP) | (sum((Tau.old - Tau)^2)/sum(Tau^2) < eps.FP)
                       
                       ## MaJ des nu_ij
                       PI       <- log(pi) - log(1-pi)
                       X2[miss] <- (logistic( log(1-psi[2]) - log(1-psi[1]) + Tau %*% PI %*% t(Tau)))[miss]
                       X3       <- 1-X2; diag(X3) <- 0
                     }
                     
                     ## Storing the current estimates
                     theta[[i]]   <- list(pi=pi, alpha=alpha)
                     psi.all[[i]] <- psi
                     Tau.all[[i]] <- Tau
                     
                     J[i] <- logLik.SBM(X2, X3, Tau, alpha, pi) - sum(Tau*log(Tau + 1*(Tau==0)))
                     if (i > 1) {
                       ## sqrt((sum((J[i] - J[i-1])^2)/sum(J[i-1]^2)))
                       conv[i] <- frobenius(theta[[i]]$pi-theta[[i-1]]$pi)/frobenius(theta[[i-1]]$pi)
                       cond    <- (i > maxIter) |  (conv[i] < eps)
                     }
                   }
                   
                 
               }
)

