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
            nBlocks          = NULL, # number of blocks 
            SBM              = NULL, # 
            sampling         = NULL, # 
            controlVEstep    = NULL,
            controlMstep     = NULL,
            ## methods
            lowerBound       = NULL, # variational lower bound (a.k.a. J)
            vICL             = NULL, # compute the (variational) integrated complete likelihood
            blocks           = NULL,  # get the most probable blocks
            # fixPoint         = NULL,
            updateNu         = NULL,
            VEstep           = function() {
              # browser()
              if(any(!(class(self$sampling) %in% c("sampling_randomPairMAR", "sampling_randomNodesMAR", "sampling_snowball")))){
                self$blockVarParam <- self$SBM$fixPoint(self$SBM, self$blockVarParam, self$completedNetwork)
                # browser()
                self$missingVarParam <- self$SBM$updateNu(self$SBM, self$sampling, self$sampledNetwork, self$blockVarParam, self$completedNetwork)
              } else {
                self$blockVarParam <- self$SBM$fixPoint_MAR(self$SBM, self$blockVarParam, self$completedNetwork, self$sampledNetwork$samplingMatrix)
              }
            },
            Mstep            = function() {
              if(any(!(class(self$sampling) %in% c("sampling_randomPairMAR", "sampling_randomNodesMAR", "sampling_snowball")))){
                # browser()
                self$SBM <- self$SBM$maximization(self$SBM, self$completedNetwork, self$blockVarParam)
              } else {
                self$SBM <- self$SBM$maximization_MAR(self$SBM, self$completedNetwork, self$blockVarParam, self$sampledNetwork$samplingMatrix)
              }
            },
            initialize       = function(SBM, sampledNetwork, sampling, nBlocks) {
              self$sampledNetwork   <- sampledNetwork
              self$SBM              <- SBM
              self$sampling         <- sampling
              self$completedNetwork <- sampledNetwork$adjacencyMatrix
              self$nBlocks          <- nBlocks
            }
          )
  )

SBM_VEMfit$set("public", "doVEM",
               function() {

                 eps        <- 1e-5
                 mc.cores   <- 2
                 zero       <- .Machine$double.eps
                 maxIter    <- 100
                 conv    <- vector("numeric", maxIter) ; conv[1] <- NA
                 theta   <- vector("list", length = maxIter)

                 cl0 <- SpectralClustering(self$completedNetwork, self$nBlocks)

                 self$blockVarParam <- matrix(0,self$SBM$nNodes,self$nBlocks) ; self$blockVarParam[cbind(1:self$SBM$nNodes, self$nBlocks)] <- 1

                 theta[[1]] <- (t(self$blockVarParam)%*% self$completedNetwork %*%self$blockVarParam) / (t(self$blockVarParam)%*%((1-diag(self$SBM$nNodes)))%*%self$blockVarParam)
                 self$completedNetwork[self$sampledNetwork$missingDyads] <- ((self$blockVarParam) %*% theta[[1]] %*% t(self$blockVarParam))[self$sampledNetwork$missingDyads]

                 i <- 0; cond <- FALSE
                 while(!cond){
                   i <- i+1
                    
                   self$Mstep()
                   # browser()
                   self$VEstep()

                   self$completedNetwork[self$sampledNetwork$missingDyads] <- self$missingVarParam
                   self$completedNetwork.bar <- 1 - self$self$completedNetwork ; diag(self$completedNetwork.bar) <- 0

                   self$lowerBound <- c(self$lowerBound, logLik.SBM(self$self$completedNetwork, self$self$completedNetwork.bar, self$blockVarParam, self$SBM$mixtureParam, self$SBM$connectParam) -
                     sum(self$blockVarParam*log(blockVarParam + 1*(blockVarParam==0))))
                   theta[[i]]      <- self$SBM$connectParam

                   if (i > 1) {
                     conv[i] <- sqrt(sum((theta[[i]]-theta[[i-1]])^2)) / sqrt(sum((theta[[i-1]])^2))
                     cond    <- (i > maxIter) |  (conv[i] < eps)
                   }
                 }
               }
)


### Tests :

# SBM : 
mySBM <- SBM_BernoulliUndirected.fit$new(100, c(1/2, 1/2), matrix(c(.2, .05, .05, .2),2,2))

# Sampled SBM :
mySampledSBM  <- sampling_doubleStandard$new(100, c(1/4, 1/2))
Y             <- mySBM$rSBM()$adjacencyMatrix
sampAdjMatrix <- mySampledSBM$rSampling(Y)$sampAdjMatrix

# Sampled network :
sampledNetwork <- sampledNetwork$new(sampAdjMatrix, FALSE)

# VEM : 
mySBM <- SBM_BernoulliUndirected.fit$new(100, NA, NA)
nBlocks <- 2
fit <- SBM_VEMfit$new(mySBM, sampledNetwork, mySampledSBM, nBlocks)
fit$doVEM()





