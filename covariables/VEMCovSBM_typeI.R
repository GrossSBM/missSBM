### ===================================================================================
### VARIATIONAL E-M FOR THE Covariate SBM (type I) CASE
### ===================================================================================

source("~/SVN/Code/code_these/functions/func_missSBM.R")

rP <- function(X,beta){
  if(is.matrix(X)) {
    return(t(X) %*% beta)
  } else {
    # browser()
      return(apply(X,c(1,2),function(x){ x %*% beta }))
    }
}

g <- function(x){
  return(-(logistic(x) - 0.5)/(0.5*x))
}

prodT <- function(X){
  if(is.matrix(X)){
    return(t(X) %*% X)
  } else {
    return(apply(X,c(1,2),function(x){ t(x) %*% x }))
  }
}


func_missSBM.CovI <- function(Y, seq.Q, cov, cl.init = "spectral", mc.cores=1, directed = FALSE){

  n <- nrow(Y)
  N <- dim(cov)[3]

  eps        <- 1e-3
  maxIter    <- 100
  eps.FP     <- 1e-2
  maxIter.FP <- 5
  mc.cores   <- 1
  zero       <- .Machine$double.eps
  if(is.character(cl.init)) stopifnot(cl.init %in% c("CAH", "spectral"))

  R  <- 1*(!is.na(Y)); diag(R) <- 1
  Y1 <- Y; Y1[is.na(Y)] <- 0
  Y0 <- (1-Y1); diag(Y0) <- 0; Y0 <- Y0 * R

  models <- mclapply(seq.Q, function(Q) {

    cl0 <- switch(cl.init,
                  "spectral"     = SpectralClustering(Y, Q),
                  "CAH"          = graphCAH(Y, Q))

    # browser()

    Tau  <- matrix(0,n,Q); Tau[cbind(1:n, cl0)] <- 1
    beta <- rnorm(N)
    pi   <- (t(Tau)%*% Y1 %*%Tau) / (t(Tau)%*%((1-diag(n))*R)%*%Tau); pi[pi > 1-zero] <- 1-zero ; pi[pi < zero] <- zero
    gamma <- log(pi)-log(1-pi)
    ksi  <- sqrt(Tau%*%(gamma^2)%*%t(Tau) + rP(cov,beta)^2 + 2*(Tau%*%gamma%*%t(Tau))*rP(cov,beta))
    # ksi(is.nan(ksi)) <- mean(ksi, na.rm=TRUE)

    ## To check the convergence : ##
    conv    <- vector("numeric", maxIter); conv[1] <- NA

    i <- 0; cond <- FALSE
    while(!cond){
      pi.old   <- pi
      i <- i+1

      gamma <- .5*(t(Tau)%*% (R*(Y1-matrix(.5,n,n)+diag(.5,n))) %*%Tau)/((t(Tau)%*% (R*(g(ksi)*(matrix(1,n,n)+2*rP(cov,beta)))) %*%Tau))
      beta  <- solve(sum(g(ksi)*prodT(cov)*R), sum(R*(Y1-matrix(.5,n,n)+diag(.5,n))) - 2*sum(R*g(ksi)*(Tau%*%(gamma^2)%*%t(Tau))))

      browser()

      pi <- 1/(1+exp(-(Tau%*%gamma%*%t(Tau) + rP(cov,beta))))
      alpha <-  colMeans(Tau)

      cond.FP <- FALSE; iter <- 0
      while(!cond.FP){
        iter <- iter + 1

        Tau <- t(gamma %*% ((R*(Y1-matrix(.5,n,n)+diag(.5,n)))%*%Tau - gamma%*%t((g(ksi)*(matrix(1,n,n) + 2*rP(cov,beta)))%*%Tau)))

        num <- rowSums(Tau)
        Tau <- Tau/num
        Tau[is.nan(Tau)] <- 0.5
        cond.FP <- (iter > maxIter.FP) | (sum((Tau.old - Tau)^2)/sum(Tau^2) < eps.FP)
      }

      ksi <- sqrt(Tau%*%(gamma^2)%*%t(Tau) + rP(cov,beta)^2 + 2*(Tau%*%(gamma^2)%*%t(Tau))*rp(cov,beta))

      if (i > 1) {
        conv[i] <- frobenius(pi-pi.old)/frobenius(pi.old)
        cond <- (i > maxIter) |  (conv[i] < eps)
      }
    }

    ICL <- -2 * (sum(Tau*log(alpha)) + .5 * sum( Y1 *(t(Tau) %*% log(pi) %*% Tau) + Y0 * (t(Tau) %*% log(1-pi) %*% Tau))) +
      Q*(Q+1)/2*log((sum(R)-n)/2) + (Q-1)*log(n)

    return(new(Class = "SBM.fit",
               theta        = list(gamma=gamma, beta = beta, alpha = alpha),
               psi          = list(NULL),
               Q            = as.integer(Q),
               cl           = lapply(Tau, function(Tau) factor(as.numeric(apply(Tau, 1, which.max)), levels = 1:Q)),
               conv         = conv[1:i],
               ICL          = ICL

    ))
  }, mc.cores=mc.cores)

  return(new(Class = "missSBM",
             missingness = "mar",
             algorithm   = "Variational Expectation-Maximization",
             obs.rate    = sum(!is.na(Y))/length(Y),
             network     = Matrix(Y),
             ICLs        = setNames(sapply(models, function(model) model@ICL),as.character(seq.Q)),
             models      = models
  ))
}
