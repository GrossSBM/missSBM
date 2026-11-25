#' Base internal class for adjusting a binary Stochastic Block Model in the context of missSBM.
#'
#' It is not designed to be called directly by the user; see the concrete variants
#' [`SimpleSBM_fit_noCov`], [`SimpleSBM_fit_withCov`] and [`SimpleSBM_fit_MNAR`].
#'
#' @import R6
SimpleSBM_fit <-
R6::R6Class(classname = "SimpleSBM_fit",
  inherit = sbm::SimpleSBM,
  private = list(
    variant = NULL, # model variant
    R       = NULL, # the sampling matrix (sparse encoding)
    S       = NULL, # the "anti" sampling matrix (sparse encoding)
    M_step  = NULL, # pointing to the appropriate M_step function
    E_step  = NULL, # pointing to the appropriate E_step function
    vLL_complete = NULL, # pointing to the appropriate Expected complete LL function
    ## memoized `imputation` (a function of Z/theta/beta only): invalidated by every method
    ## that can change any of those, recomputed lazily on the next read
    imputation_cache = NULL
  ),
  public = list(
    #' @description constructor for simpleSBM_fit for missSBM purpose
    #' @param networkData a structure to store network under missing data condition: either a matrix possibly with NA, or a missSBM:::partlyObservedNetwork
    #' @param clusterInit Initial clustering: a vector with size \code{ncol(adjacencyMatrix)}, providing a user-defined clustering with \code{nbBlocks} levels.
    #' @param covarList An optional list with M entries (the M covariates).
    initialize = function(networkData, clusterInit, covarList = list()) {

      ## networkData
      if (inherits(networkData, "matrix")) networkData <- partlyObservedNetwork$new(networkData, covarList)
      stopifnot(inherits(networkData, "partlyObservedNetwork"))

      ## Initial Clustering
      private$Z <- clustering_indicator(clusterInit)

      # Basic fields initialization by call to super constructor
### TODO: determine the model according to the adjacency matrix
      model <- "bernoulli"
      super$initialize(model        = model,
                       directed     = networkData$is_directed,
                       nbNodes      = networkData$nbNodes,
                       blockProp    = vector("numeric", ncol(private$Z)) + .Machine$double.eps,
### FIXME: this should be model-dependent
                       connectParam = list(mean = matrix(0, ncol(private$Z), ncol(private$Z))),
                       covarList    = covarList)

      # Storing data
      private$R <- networkData$samplingMatrix
      private$S <- networkData$samplingMatrixBar
      private$Y <- networkData$networkData

      ## point to the functions that performs E/M steps and compute the likelihood
      directedSuffix <- ifelse(self$directed, "directed", "undirected")
      covarSuffix    <- ifelse(self$nbCovariates > 0, "covariates", "nocovariate")
      private$variant       <- paste(model, directedSuffix, covarSuffix, sep = "_")
      private$M_step        <- get(paste("M_step_sparse",       model, covarSuffix, sep = "_"))
      private$E_step        <- get(paste("E_step_sparse",       model, covarSuffix, sep = "_"))
      private$vLL_complete  <- get(paste("vLL_complete_sparse", model, covarSuffix, sep = "_"))

### TODO:
###  - check if parameters are not already initialize
###  - specialize the initialization to each model (this should be done in sbm::SimpleSBM...)

      ## Initialize estimation of the model parameters
      private$theta$mean <- matrix(.Machine$double.eps, ncol(private$Z), ncol(private$Z))
      private$beta       <- numeric(self$nbCovariates)
      self$update_parameters()

      invisible(self)
    },
    #' @description method to perform estimation via variational EM
    #' @param threshold stop when an optimization step changes the objective function by less than threshold. Default is 1e-4.
    #' @param maxIter V-EM algorithm stops when the number of iteration exceeds maxIter. Default is 10
    #' @param fixPointIter number of fix-point iterations in the Variational E step. Default is 5.
    #' @param trace logical for verbosity. Default is \code{FALSE}.
    doVEM = function(threshold = 1e-2, maxIter = 100, fixPointIter = 3, trace = FALSE) {

      if (trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")

      run_VEM(
        control        = list(threshold = threshold, maxIter = maxIter, fixPointIter = fixPointIter, trace = trace),
        init_stop      = self$nbBlocks <= 1,
        e_step         = function(fixPointIter) for (i in seq.int(fixPointIter)) self$update_blocks(),
        m_step         = function() self$update_parameters(),
        get_loglik     = function() self$loglik,
        get_theta      = function() private$theta$mean,
        snapshot       = function() self$get_state(),
        restore        = function(state) self$set_state(state),
        reorder        = function() self$reorder()
      )
    },
    #' @description permute group labels by order of decreasing probability
    #' @return invisibly, the permutation \code{o} applied (new block \code{k} is former block
    #'   \code{o[k]}), so that the parameters of the sampling indexed by the blocks can follow
    reorder = function(){
      private$imputation_cache <- NULL
      o <- order(private$theta$mean %*% private$pi, decreasing = TRUE)
      private$pi <- private$pi[o]
      private$theta$mean <- private$theta$mean[o, o, drop = FALSE]
      private$Z <- private$Z[, o, drop = FALSE]
      invisible(o)
    },
    #' @description a lightweight snapshot of the mutable VEM state (as opposed to \code{clone()},
    #' which duplicates the whole object, including the -- possibly large -- network data)
    get_state = function() {
      list(theta = private$theta, Z = private$Z, pi = private$pi, beta = private$beta)
    },
    #' @description restore a state previously returned by \code{get_state()}
    #' @param state a state, as returned by \code{get_state()}
    set_state = function(state) {
      private$imputation_cache <- NULL
      private$theta <- state$theta
      private$Z     <- state$Z
      private$pi    <- state$pi
      private$beta  <- state$beta
      invisible(self)
    }
  ),
  active = list(
    #' @field type the type of SBM (distribution of edges values, network type, presence of covariates)
    type = function(value) {private$variant},
    #' @field nbObservedDyads the number of observed dyads
    nbObservedDyads = function(value) {sum(private$R)},
    #' @field penalty double, value of the penalty term in ICL: the connectivity and the effects
    #'   of the covariates are estimated on the observed dyads, which are counted instead of all
    #'   the possible dyads (as sbm does for a network with missing values)
    penalty  = function(value) {unname((self$nbConnectParam + self$nbCovariates) * log(self$nbObservedDyads) + (self$nbBlocks-1) * log(self$nbNodes))},
    #' @field entropy double, value of the entropy due to the clustering distribution
    entropy  = function(value) {-sum(xlogx(private$Z))},
    #' @field loglik double: approximation of the log-likelihood (variational lower bound) reached
    loglik = function(value) {self$vExpec + self$entropy},
    #' @field ICL double: value of the integrated classification log-likelihood: -2 times the
    #'   variational expectation of the complete log-likelihood at the MAP memberships (each node
    #'   in its most probable block), plus the penalty, as for [missSBM_fit]
    ICL    = function(value) {
      tau <- private$Z
      Z <- matrix(0, nrow(tau), ncol(tau))
      Z[cbind(seq_len(nrow(tau)), max.col(tau, ties.method = "first"))] <- 1
      -2 * self$vExpec_at(Z) + self$penalty
    }
  )
)

#' Internal class for a binary SBM fit on the observed dyads, without covariates: MAR sampling designs, 'block-node' and 'block-dyad' (with their sampling term in the update of the blocks) and 'degree'.
#'
#' It is not designed to be called directly by the user.
#'
#' @import R6
SimpleSBM_fit_noCov <-
R6::R6Class(classname = "SimpleSBM_fit_noCov",
  inherit = SimpleSBM_fit,
  public = list(
    #' @description update parameters estimation (M-step)
    #' @param ... additional arguments, only required for MNAR cases
    update_parameters = function(...) {
      private$imputation_cache <- NULL
      res <- private$M_step(private$Y, private$R, private$Z, !self$directed)
      private$theta <- res$theta
      private$theta$mean <- check_boundaries(private$theta$mean)
      private$pi    <- check_boundaries(as.numeric(res$pi))
      invisible(res)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda term of the sampling design depending on the blocks ('block-node' and
    #'   'block-dyad'), added to the update of the blocks; 0 for the other designs
    #' @param ... additional arguments, only required for MNAR cases
    update_blocks =   function(log_lambda = 0, ...) {
      private$imputation_cache <- NULL
      if (identical(log_lambda, 0)) {
        private$Z <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi)
      } else if (self$nbBlocks > 1) {
        ## sampling depending on the blocks (block-node, block-dyad): observed dyads + sampling term
        log_tau   <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi, rescale = FALSE)
        private$Z <- t(apply(log_tau + log_lambda, 1, .softmax))
      }
    },
    #' @description for each node, the complete-data log-likelihood it would contribute to
    #'   each class if hard-assigned there (theta/pi held fixed), used to decide node-swap
    #'   moves in \code{missSBM_fit$polish()}.
    #' @param log_lambda additional sampling-design-dependent term, added as-is
    #' @param ... additional arguments, not used
    #' @return an N x Q matrix
    polish_log_tau = function(log_lambda = 0, ...) {
      private$E_step(private$Y, private$R, self$indMemberships, private$theta$mean, private$pi, rescale = FALSE) + log_lambda
    },
    #' @description variational expectation of the log-likelihood of the missing dyads, given
    #'   their imputations (used by the ICL of 'degree', see [missSBM_fit])
    #' @param nu the imputations of the missing dyads (sparse matrix with the pattern of the
    #'   missing dyads)
    #' @param Z an n x Q matrix of memberships (default: the current variational ones)
    vExpec_imputed = function(nu, Z = private$Z) {
      private$vLL_complete(nu, private$S, Z, private$theta$mean, private$pi) - sum(Z %*% log(private$pi))
    },
    #' @description variational expectation of the log-likelihood of the observed dyads and of
    #'   the blocks at given memberships (used by the ICL, at the MAP memberships)
    #' @param Z an n x Q matrix of memberships
    vExpec_at = function(Z) {
      private$vLL_complete(private$Y, private$R, Z, private$theta$mean, private$pi)
    }
  ),
  active = list(
    #' @field imputation the matrix of imputed values
    imputation = function(value) {
      if (is.null(private$imputation_cache))
        private$imputation_cache <- .mask_dense_at_pattern(private$Z %*% log(private$theta$mean/(1-private$theta$mean)) %*% t(private$Z), private$S, .logistic)
      private$imputation_cache
    },
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      private$vLL_complete(private$Y, private$R, private$Z, private$theta$mean, private$pi)
    }
  )
)

#' Internal class for a binary SBM fit under MAR sampling designs with covariates.
#'
#' It is not designed to be called directly by the user.
#'
#' @import R6
SimpleSBM_fit_withCov <-
R6::R6Class(classname = "SimpleSBM_fit_withCov",
  inherit = SimpleSBM_fit,
  public = list(
    #' @description update parameters estimation (M-step) via Newton-Raphson: the M-step
    #'   objective is a weighted logistic regression (concave), so Newton converges in a
    #'   handful of iterations -- no external optimizer is required.
    #' @param ... use for compatibility
    update_parameters = function(...) {
      private$imputation_cache <- NULL
      res <- private$M_step(
        init_param = list(Gamma = .logit(private$theta$mean), beta = private$beta),
        Y = private$Y,
        R = private$R,
        X = self$covarArray,
        Z = private$Z,
        !self$directed,
        maxIter = 50,
        tol = 1e-10
      )
      private$beta  <- as.numeric(res$beta)
      private$theta <- res$theta
      private$pi    <- as.numeric(res$pi)
      invisible(res)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param ... use for compatibility
    update_blocks =   function(...) {
      private$imputation_cache <- NULL
      private$Z <- private$E_step(private$Y, private$R, self$covarEffect, private$Z, .logit(private$theta$mean), private$pi, !self$directed, TRUE)
    },
    #' @description for each node, the complete-data log-likelihood it would contribute to
    #'   each class if hard-assigned there (theta/beta/pi held fixed), used to decide node-swap
    #'   moves in \code{missSBM_fit$polish()}.
    #' @param log_lambda additional sampling-design-dependent term, added as-is
    #' @param ... additional arguments, not used
    #' @return an N x Q matrix
    polish_log_tau = function(log_lambda = 0, ...) {
      private$E_step(private$Y, private$R, self$covarEffect, self$indMemberships, .logit(private$theta$mean), private$pi, !self$directed, FALSE) + log_lambda
    },
    #' @description variational expectation of the log-likelihood of the observed dyads and of
    #'   the blocks at given memberships (used by the ICL, at the MAP memberships)
    #' @param Z an n x Q matrix of memberships
    vExpec_at = function(Z) {
      private$vLL_complete(private$Y, private$R, self$covarEffect, Z, .logit(private$theta$mean), private$pi)
    }
  ),
  active = list(
    #' @field imputation the matrix of imputed values
    imputation = function(value) {
      if (is.null(private$imputation_cache))
        private$imputation_cache <- .mask_dense_at_pattern(private$Z %*% .logit(private$theta$mean) %*% t(private$Z) + self$covarEffect, private$S, .logistic)
      private$imputation_cache
    },
    #' @field vExpec double: variational approximation of the expectation complete log-likelihood
    vExpec = function(value) {
      private$vLL_complete(private$Y, private$R, self$covarEffect, private$Z, .logit(private$theta$mean), private$pi)
    }
  )
)

#' Internal class for a binary SBM fit under the double-standard sampling design (MNAR,
#' depending on the value of each dyad).
#'
#' The missing dyads are integrated out given the blocks. A dyad between blocks q and l is
#' missing with probability \eqn{m_{ql} = \pi_{ql}(1 - p_1) + (1 - \pi_{ql})(1 - p_0)}, where
#' \eqn{(p_0, p_1)} are the parameters of the sampling: the update of the blocks gets the term
#' \eqn{\sum_{j} \sum_l \tau_{jl} \log m_{ql}} over the missing dyads (i, j). The missing dyads
#' tie the connectivity to \eqn{(p_0, p_1)}: the M-step maximizes over the three jointly, by EM
#' on the missing dyads (the posterior probability that a missing dyad between blocks q and l
#' is an edge is \eqn{\nu_{ql} = \pi_{ql}(1 - p_1) / m_{ql}}), iterated to convergence, all in
#' Q x Q given the memberships.
#'
#' It is not designed to be called directly by the user.
#'
#' @import R6
SimpleSBM_fit_MNAR <-
R6::R6Class(classname = "SimpleSBM_fit_MNAR",
  inherit = SimpleSBM_fit_noCov,
  private = list(
    ## term of the missing dyads in the update of the blocks: for each node i and block q, the sum
    ## over its missing dyads (i, j) of sum_l Z_jl log m_ql (both triangles if undirected)
    missing_term = function(Z, log_m) {
      if (identical(log_m, 0)) return(0) # no sampling term (the SBM alone)
      as.matrix(private$S %*% Z %*% t(log_m) + Matrix::crossprod(private$S, Z) %*% log_m) # no transpose of S (#49)
    }
  ),
  public = list(
    #' @description update parameters estimation (M-step)
    #' @param psi the parameters (p0, p1) of the sampling, starting point of the joint M-step.
    #'   If \code{NULL}, or with a single block (where the sampling is constrained to dyad
    #'   sampling), the SBM is fitted on the observed dyads.
    #' @param tol,maxit tolerance on the parameters and maximal number of iterations of the EM
    #'   on the missing dyads
    #' @return invisibly, the updated (p0, p1), or \code{NULL} if the SBM is fitted on the
    #'   observed dyads
    update_parameters = function(psi = NULL, tol = 1e-8, maxit = 10000) {
      private$imputation_cache <- NULL
      if (is.null(psi) || self$nbBlocks == 1) {
        super$update_parameters()
        return(invisible(NULL))
      }
      Z <- private$Z
      Y1 <- as.matrix(t(Z) %*% private$Y %*% Z) # observed edges between blocks q and l
      N  <- as.matrix(t(Z) %*% private$R %*% Z) # observed dyads
      M  <- as.matrix(t(Z) %*% private$S %*% Z) # missing dyads
      ## the statistics of the sampling count each dyad once
      So <- sum(Y1); So.bar <- sum(N) - So; Dm <- sum(M)
      if (!self$directed) {Y1 <- Y1 + t(Y1); N <- N + t(N); M <- M + t(M)}
      theta <- private$theta$mean
      for (k in seq_len(maxit)) {
        theta_old <- theta; psi_old <- psi
        ## E: posterior probability that a missing dyad between blocks q and l is an edge
        nu <- theta * (1 - psi[2]) / (theta * (1 - psi[2]) + (1 - theta) * (1 - psi[1]))
        ## M: the missing dyads completed by nu
        theta <- check_boundaries((Y1 + nu * M) / (N + M))
        Sm  <- Dm * sum(nu * M) / sum(M) # expected edges among the missing dyads, each counted once
        psi <- check_boundaries(c(So.bar / (So.bar + Dm - Sm), So / (So + Sm)))
        if (max(abs(theta - theta_old), abs(psi - psi_old)) < tol) break
      }
      private$theta$mean <- theta
      private$pi <- check_boundaries(colMeans(Z))
      invisible(psi)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda Q x Q matrix, the log-probability that a dyad between blocks q and l is
    #'   missing (field \code{log_lambda} of [doubleStandardSampling_fit])
    update_blocks = function(log_lambda = 0) {
      private$imputation_cache <- NULL
      if (self$nbBlocks > 1) {
        log_tau <- private$E_step(private$Y, private$R, private$Z, private$theta$mean, private$pi, rescale = FALSE) +
          private$missing_term(private$Z, log_lambda)
        private$Z <- t(apply(log_tau, 1, .softmax))
      }
    },
    #' @description for each node, the log-likelihood of the observed data it would contribute
    #'   to each class if hard-assigned there (parameters held fixed), used to decide node-swap
    #'   moves in \code{missSBM_fit$polish()}.
    #' @param log_lambda Q x Q matrix, the log-probability that a dyad between blocks q and l is
    #'   missing
    #' @param ... additional arguments, not used
    #' @return an N x Q matrix
    polish_log_tau = function(log_lambda = 0, ...) {
      Zhard <- self$indMemberships
      private$E_step(private$Y, private$R, Zhard, private$theta$mean, private$pi, rescale = FALSE) +
        private$missing_term(Zhard, log_lambda)
    }
  )
)

#' Internal class for a binary SBM fit under the degree sampling design (MNAR, depending on the
#' degrees).
#'
#' Without imputations (\code{nu = NULL}), the SBM is fitted on the observed dyads, as
#' [`SimpleSBM_fit_noCov`]: the first stage of the estimation of this design. With the
#' imputations of the sampling, the blocks and the connectivity are updated on the network
#' completed by them: the coordinate ascent of the variational bound of the degree sampling
#' (factorized q(Z) q(Y^m)), its third stage (see [missSBM_fit]). Only the pattern of the missing
#' dyads is used for the imputations.
#'
#' It is not designed to be called directly by the user.
#'
#' @import R6
SimpleSBM_fit_degree <-
R6::R6Class(classname = "SimpleSBM_fit_degree",
  inherit = SimpleSBM_fit_noCov,
  private = list(
    ## log tau (up to a constant) on the network completed by nu: the E-step is linear in (Y, R),
    ## so it is the sum of the E-steps of the observed and of the imputed dyads, log(pi) counted once
    log_tau_completed = function(Z, nu) {
      private$E_step(private$Y, private$R, Z, private$theta$mean, private$pi, rescale = FALSE) +
        private$E_step(nu, private$S, Z, private$theta$mean, private$pi, rescale = FALSE) -
        matrix(log(private$pi), nrow(Z), ncol(Z), byrow = TRUE)
    }
  ),
  public = list(
    #' @description update parameters estimation (M-step)
    #' @param nu the imputations of the missing dyads (sparse matrix with the pattern of the
    #'   missing dyads), or \code{NULL} to fit the SBM on the observed dyads
    #' @param ... additional arguments, not used
    update_parameters = function(nu = NULL, ...) {
      if (is.null(nu)) return(super$update_parameters())
      private$imputation_cache <- NULL
      Z  <- private$Z
      N1 <- as.matrix(t(Z) %*% private$Y %*% Z + t(Z) %*% nu %*% Z)
      N  <- as.matrix(t(Z) %*% private$R %*% Z + t(Z) %*% private$S %*% Z)
      if (!self$directed) {N1 <- N1 + t(N1); N <- N + t(N)}
      private$theta$mean <- check_boundaries(N1 / N)
      private$pi <- check_boundaries(colMeans(Z))
      invisible(self)
    },
    #' @description update variational estimation of blocks (VE-step)
    #' @param log_lambda not used (0 for the degree sampling)
    #' @param nu the imputations of the missing dyads, or \code{NULL} to use the observed dyads only
    #' @param ... additional arguments, not used
    update_blocks = function(log_lambda = 0, nu = NULL, ...) {
      if (is.null(nu)) return(super$update_blocks(log_lambda))
      private$imputation_cache <- NULL
      if (self$nbBlocks > 1) private$Z <- t(apply(private$log_tau_completed(private$Z, nu), 1, .softmax))
    },
    #' @description for each node, the expected complete log-likelihood it would contribute to
    #'   each class if hard-assigned there (parameters and imputations held fixed), used to decide
    #'   node-swap moves in \code{missSBM_fit$polish()}.
    #' @param log_lambda not used (0 for the degree sampling)
    #' @param nu the imputations of the missing dyads, or \code{NULL} to use the observed dyads only
    #' @param ... additional arguments, not used
    #' @return an N x Q matrix
    polish_log_tau = function(log_lambda = 0, nu = NULL, ...) {
      if (is.null(nu)) return(super$polish_log_tau(log_lambda))
      private$log_tau_completed(self$indMemberships, nu)
    }
  )
)
