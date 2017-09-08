#' a collection of adjusted Stochastic Block Model
#'
#' @importFrom R6 R6Class
SBM_collection <-
  R6Class(classname = "SBM_collection",
    public = list(
      sampledNetwork = NULL, # the sampled network data
      vBlocks  = NULL, # the vector of number of blocks considered in the collection
      models   = NULL, # the collection of SBM fit
      sampling = NULL, # the sampling design for missing data modeling
      samplingData = NULL, # informations about the sampling
      family   = NULL, # the emission law of the adjacency matrix
      link     = NULL  # directed or not, depends on what we are modeling
      ## methods
      ### TODO: VE-step and M-step
      ## criteria = NULL, # variational bound and ICL for all models
  )
)

SBM_collection$set("public", "initialize",
  function(sampledNetwork, vBlocks, sampling, family, link) {
    self$family   <- family
    self$link     <- link
    self$vBlocks  <- vBlocks
    self$sampledNetwork <- sampledNetwork$new(sampledNetwork, link)
    self$samplingData       <- switch(sampling,
                            "doubleStandard" = sampling_doubleStandard$new(self$sampledNetwork$nNodes, NA, link),
                            "class"          = sampling_class$new(self$sampledNetwork$nNodes, NA, link),
                            "starDegree"     = sampling_starDegree$new(self$sampledNetwork$nNodes, NA, link),
                            "MAREdge"        = sampling_randomPairMAR$new(self$sampledNetwork$nNodes, NA, link),
                            "MARNode"        = sampling_randomNodeMAR$new(self$sampledNetwork$nNodes, NA, link),
                            "snowball"       = sampling_snowball$new(self$sampledNetwork$nNodes, NA, link))
    self$sampling <- sampling
  
  }
)

SBM_collection$set("public", "estimate",
  function() {

    ## A. INITIALISATION DES PARAMÈTRES
    
    ## A.1 PARAMÈTRES CONTRÔLANT L'OPTIMISATION
    eps        <- 1e-5
    maxIter    <- 250
    eps.FP     <- 1e-2
    maxIter.FP <- 5
    mc.cores   <- mc.cores
    zero       <- .Machine$double.eps

    ## A.2 PARAMÈTRES FIXES EN COURS D'ALGORITHME (données)

    SBM <- switch(paste0("SBM_",self$family, self$link),
                  "SBM_bernoulliUndirected" = SBM_bernoulliUndirected$new(self$sampling$nNodes, NA, NA),
                  "SBM_bernoulliDirected"   = SBM_bernoulliDirected$new(self$sampling$nNodes, NA, NA),
                  "SBM_poissonUndirected"   = SBM_poissonUndirected$new(self$sampling$nNodes, NA, NA),
                  "SBM_poissonDirected"     = SBM_poissonDirected$new(self$sampling$nNodes, NA, NA))
    

    VEMfit <- SBM_VEMfit$new(self$sampledNetwork$adjacencyMatrix, NA, NA, SBM, lowerBound, vICL, blocks)

    
    
    self$models <- mclapply(self$vBlocks, function(Q) {
      ## ===================================================================
      ## B. VEM À NOMBRE DE CLASSES Q FIXÉ
      
      ## B.1 INITIALISATION DES CLUSTERS
      # if(is.null(CL)){
      #   if (is.character(cl.init) | is.null(cl.init)) { ## CAH, spectral clustering or user defined
      #     cl0 <- switch(cl.init,
      #                   "spectral"     = SpectralClustering(X, Q),
      #                   "CAH"          = graphCAH(X, Q),
      #                   "mar.spectral" = SpectralClustering(X, Q),
      #                   "mar.CAH"      = graphCAH(X, Q))
      #   } else {
      #     if (is.list(cl.init)) {
      #       cl0 <- cl.init[[Q]]
      #     } else {
      #       cl0 <- cl.init[[which(seq.Q == Q)]]
      #     }
      #   }
      # } else {
      #   cl0 <- CL
      # }
      # if(!is.null(CL.init) & is.null(CL)){
      #   cl0 <- CL.init
      # }
      
      
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
      
      ## Calcul de l'ICL :
      
      logLikFull <- logLik.SBM(X2, X3, Tau, alpha, pi) + logLik.missing.doublestandard.VEM(X2, miss, obs, psi, Tau)
      
      ICL <- -2*logLikFull + (2 + Q*(Q+1)/2)*log(n*(n-1)/2) + (Q-1)*log(n)
      
      return()
    }, mc.cores=mc.cores)
    
    return()
    
  }
)

SBM_collection$set("public", "doVEstep",
  function(model) {

  }
)

SBM_collection$set("public", "doMstep",
  function(model) {

  }
)
