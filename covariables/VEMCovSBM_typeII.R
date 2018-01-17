### ===================================================================================
### VARIATIONAL E-M FOR THE Covariate SBM (type II) CASE
### ===================================================================================

source("~/SVN/Code/code_these/functions/func_missSBM.R")

func_missSBM.CovII <- function(Y, seq.Q, cov, cl.init = "spectral", mc.cores=1, directed = FALSE){

  n <- nrow(Y)
  N <- nrow(cov)

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

    Tau  <- matrix(0,n,Q); Tau[cbind(1:n, cl0)] <- 1
    beta <- matrix(rnorm(N*(Q-1)), nrow=N, ncol=Q-1)

    theta   <- vector("list", length = maxIter)
    Tau.all <- vector("list", length = maxIter)
    conv    <- vector("numeric", maxIter); conv[1] <- NA

    i <- 0; cond <- FALSE
    while(!cond){
      i <- i+1
      pi <- (t(Tau)%*% Y1 %*%Tau) / (t(Tau)%*%((1-diag(n))*R)%*%Tau)
      pi[is.nan(pi)] <- zero ; pi[pi > 1-zero] <- 1-zero ; pi[pi < zero] <- zero

      ### MAJ alpha/beta :

      if(Q != 1){
        fr <- function(x) {
          b <- matrix(x, nrow=N, ncol=Q-1, byrow=F)
          a <- drawAlpha(n,Q,b,cov)
          return(sum(Tau*log(a)))
        }
        grr <- function(x) {
          b <- matrix(x, nrow=N, ncol=Q-1, byrow=F)
          grad <- NULL
          for(q in 1:(Q-1)){
            acc <- rep(0,N)
            for(i in 1:n){
              acc <- acc + cov[,i]*Tau[i,q]*(1 - exp(t(beta[,q]) %*% cov[,i])/(1 + sum(exp(t(beta) %*% cov[,i]))))
            }
            grad <- cbind(grad, acc)
          }
          return(as.numeric(grad))
        }
        Optim <- optim(beta, fr, grr, method = "BFGS", control = list(fnscale = -1))

        beta <- matrix(Optim$par, nrow=N, ncol=Q-1, byrow=F)
        alpha <-  drawAlpha(n,Q,beta,cov)
      } else {
        alpha <-  colMeans(Tau)
      }

      ###

      cond.FP <- FALSE; iter <- 0
      while(!cond.FP){
        Tau.old <- Tau
        iter <- iter + 1
        if(!directed){
          Tau <- exp(sweep(Y1 %*% Tau %*% t(log(pi)) + Y0 %*% Tau %*% t(log(1-pi)) + t(Y1) %*% Tau %*% log(pi) + t(Y0) %*% Tau %*% log(1-pi),2,log(alpha),"+"))
        } else {
          Tau <- exp(sweep(Y1 %*% Tau %*% t(log(pi)) + Y0 %*% Tau %*% t(log(1-pi)),2,log(alpha),"+"))
        }
        num <- rowSums(Tau)
        Tau <- Tau/num
        Tau[is.nan(Tau)] <- 0.5
        cond.FP <- (iter > maxIter.FP) | (sum((Tau.old - Tau)^2)/sum(Tau^2) < eps.FP)
      }

      theta[[i]]   <- list(pi=pi, beta = beta, alpha = alpha)
      Tau.all[[i]] <- Tau

      if (i > 1) {
        conv[i] <- frobenius(theta[[i]]$pi-theta[[i-1]]$pi)/frobenius(theta[[i-1]]$pi)
        cond <- (i > maxIter) |  (conv[i] < eps)
      }
    }

    ICL <- -2 * (sum(Tau*log(alpha)) + .5 * sum( Y1 *(Tau %*% log(pi) %*% t(Tau)) + Y0 * (Tau %*% log(1-pi) %*% t(Tau)))) +
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
