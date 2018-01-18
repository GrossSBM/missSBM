### ===================================================================================
### VARIATIONAL E-M FOR THE Covariate SBM (type I) CASE
### ===================================================================================

source("~/SVN/Code/code_these/functions/func_missSBM.R")

rP <- function(X,beta){
  if(is.matrix(X)) {
    return(t(X) %*% beta)
  } else {
      return(apply(X,c(1,2),function(x){ x %*% beta }))
    }
}

g <- function(x){
  return(-(logistic(x) - 0.5)/(0.5*x))
}

# prodT <- function(X){
#   return(apply(X,c(1,2),function(x){ x %*% t(x) }))
# }

sysA <- function(X,ksi,R){
  K <- dim(X)[3]
  acc <- matrix(0,K,K)
  for (i in 1:n) {
    for (j in 1:n) {
      acc <- acc + (X[i,j,] %*% t(X[i,j,]))*ksi[i,j]*R[i,j]
    }
  }
  return(acc)
}

sysB <- function(Y,ksi,Tau,gamma,X,R){
  n <- dim(X)[1]
  K <- dim(X)[3]
  acc <- rep(0,K)
  for (i in 1:n) {
    for (j in 1:n) {
      acc <- acc + R[i,j]*(Y[i,j] - .5 - 2*g(ksi[i,j])*(Tau[i,]%*%gamma%*%Tau[j,]))*X[i,j,]
    }
  }
  return(acc)
}


func_missSBM.CovI <- function(Y, seq.Q, cov, cl.init = "spectral", mc.cores=1, directed = FALSE){

  n <- nrow(Y)
  N <- dim(cov)[3]

  eps        <- .4
  maxIter    <- 50
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

    Tau  <- matrix(0,n,Q); Tau[cbind(1:n, cl0)] <- 1
    beta <- rnorm(N)
    pi   <- (t(Tau)%*% Y1 %*%Tau) / (t(Tau)%*%((1-diag(n))*R)%*%Tau); pi[pi > 1-zero] <- 1-zero ; pi[pi < zero] <- zero
    gamma <- log(pi)-log(1-pi)
    ksi  <- sqrt((Tau%*%gamma%*%t(Tau) + rP(cov,beta))^2)

    Tau.all <- vector("list", length = maxIter)
    theta   <- vector("list", length = maxIter)


    ## To check the convergence : ##
    conv    <- vector("numeric", maxIter); conv[1] <- NA

    i <- 0; cond <- FALSE
    while(!cond){
      pi.old   <- pi
      i <- i+1

      gamma <- .5*(t(Tau)%*% (R*(Y1-matrix(.5,n,n)+diag(.5,n))) %*%Tau)/((t(Tau)%*% (R*(g(ksi)*(matrix(1,n,n)+2*rP(cov,beta)))) %*%Tau))
      beta  <- solve(sysA(cov,ksi,R), sysB(Y1,ksi,Tau,gamma,cov,R))

      pi <- 1/(1+exp(-(Tau%*%gamma%*%t(Tau) + rP(cov,beta))))
      alpha <-  colMeans(Tau)

      cond.FP <- FALSE; iter <- 0
      while(!cond.FP){
        iter <- iter + 1
        Tau.old <- Tau

        # browser()
        Tau <-exp(sweep((Y1-matrix(.5,n,n)+diag(.5,n))%*%Tau %*%t(gamma) - (R*g(ksi)*(matrix(1,n,n) + 2*rP(cov,beta)))%*%Tau%*%t(gamma^2) ,2,log(alpha),"+"))

        num <- rowSums(Tau)
        Tau <- Tau/num
        Tau[is.nan(Tau)] <- 0.5
        cond.FP <- (iter > maxIter.FP) | (sum((Tau.old - Tau)^2)/sum(Tau^2) < eps.FP)
      }

      # browser()
      ksi <- sqrt((Tau%*%gamma%*%t(Tau) + rP(cov,beta))^2)
      Tau.all[[i]] <- Tau
      theta[[i]]   <- list(pi=pi, alpha=alpha, gamma=gamma)


      if (i > 1) {
        conv[i] <- frobenius(pi-pi.old)/frobenius(pi.old)
        print(conv[i])
        cond <- (i > maxIter) |  (conv[i] < eps)
      }
    }
    # browser()

    ICL <- -2 * (sum(Tau*log(alpha)) + .5 * sum( Y1*log(pi) + Y0*log(1-pi))) +
      Q*(Q+1)/2*log((sum(R)-n)/2) + (Q-1)*log(n)

    return(new(Class = "SBM.fit",
               theta        = theta[i],
               psi          = list(NULL),
               Q            = as.integer(Q),
               cl           = lapply(Tau.all[i], function(Tau) factor(as.numeric(apply(Tau, 1, which.max)), levels = 1:Q)),
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
