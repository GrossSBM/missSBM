#' Virtual class used to define a family of networkSamplingDyads_fit
#' @include R6Class-networkSampling.R
#' @import R6
networkSamplingDyads_fit <-
  R6::R6Class(classname = "networkSamplingDyads_fit",
  inherit = networkSampling,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    card_D   = NULL, # number of possible dyads in the network
    card_D_o = NULL, # number of observed dyads in the network
    card_D_m = NULL  # number of missing dyads in the network
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description constructor for networkSampling_fit
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param name a character for the name of sampling to fit on the partlyObservedNetwork
    initialize = function(partlyObservedNetwork, name) {
      private$name     <- name
      private$card_D   <- partlyObservedNetwork$nbDyads
      private$card_D_o <- sum(partlyObservedNetwork$samplingMatrix)
      private$card_D_m <- private$card_D - private$card_D_o
    },
    #' @description show method
    show = function() {
      super$show()
      cat("  $penalty, $vExpec\n")
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param ... use for compatibility
    update_parameters = function(...) {invisible(NULL)},
    #' @description a method to update the imputation of the missing entries.
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    update_imputation = function(nu) {nu}, # good for MCAR on node, dyads and MNAR with blocks
    #' @description permute the parameters indexed by the blocks, after the blocks of the SBM
    #'   were permuted. By default, nothing to do (no parameter indexed by the blocks)
    #' @param o the permutation applied to the blocks (new block k is former block \code{o[k]})
    reorder = function(o) {invisible(self)}
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## ACTIVE BINDING
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field penalty double, value of the penalty term in ICL
    penalty = function(value) {log(private$card_D) * self$df},
    #' @field log_lambda double, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {0}
  )
)

#' Virtual class used to define a family of networkSamplingNodes_fit
#' @include R6Class-networkSampling.R
#' @import R6
networkSamplingNodes_fit <-
  R6::R6Class(classname = "networkSamplingNodes_fit",
  inherit = networkSampling,
  private = list(
    card_N   = NULL, # number of nodes in the network
    card_N_m = NULL, # number of observed nodes
    card_N_o = NULL, # number of missing nodes
    N_obs    = NULL  # boolean for observed nodes
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param name a character for the name of sampling to fit on the partlyObservedNetwork
    initialize = function(partlyObservedNetwork, name) {
      private$name     <- name
      private$N_obs    <- partlyObservedNetwork$observedNodes
      private$card_N   <- partlyObservedNetwork$nbNodes
      private$card_N_o <- sum(partlyObservedNetwork$observedNodes)
      private$card_N_m <- private$card_N - private$card_N_o
    },
    #' @description show method
    show = function() {
      super$show()
      cat("  $penalty, $vExpec\n")
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param ... use for compatibility
    update_parameters = function(...) {invisible(NULL)},
    #' @description a method to update the imputation of the missing entries.
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    update_imputation = function(nu) {nu}, ## good for MCAR on node, dyads and MNAR with blocks
    #' @description permute the parameters indexed by the blocks, after the blocks of the SBM
    #'   were permuted. By default, nothing to do (no parameter indexed by the blocks)
    #' @param o the permutation applied to the blocks (new block k is former block \code{o[k]})
    reorder = function(o) {invisible(self)}
  ),
  active = list(
    #' @field penalty double, value of the penalty term in ICL
    penalty = function(value) {log(private$card_N) * self$df},
    #' @field log_lambda double, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {0}
  )
)

#' Class for fitting a dyad sampling
dyadSampling_fit <-
  R6::R6Class(classname = "dyadSampling_fit",
  inherit = networkSamplingDyads_fit,
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "dyad")
      private$psi <- check_boundaries(private$card_D_o / private$card_D)
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      res <- private$card_D_o * log(private$psi + 1 * (private$psi == 0)) +
        private$card_D_m * log(1 - private$psi + 1 * (private$psi == 1))
      res
    }
  )
)

#' Class for fitting a dyad sampling with covariates
covarDyadSampling_fit <-
  R6::R6Class(classname = "covarDyadSampling_fit",
  inherit = networkSamplingDyads_fit,
  public = list(
    #' @description constructor
    #' @param partialNet a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(partialNet, ...) {
      super$initialize(partialNet, "covar-dyad")
      dyads <- which(valid_dyads(dim(partialNet$networkData), partialNet$is_directed), arr.ind = TRUE)
      X <- cbind(1, apply(partialNet$covarArray, 3, function(x) x[dyads]))
      y <- partialNet$samplingMatrix[dyads]
      glm_out     <- suppressWarnings(glm.fit(X, y, family = binomial()))
      private$psi <- coefficients(glm_out)
      y_hat <- fitted(glm_out)
      private$rho <- list(obs = y_hat[y == 1], miss = y_hat[y == 0])
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      res <- sum(log(1 - private$rho$miss)) + sum(log(private$rho$obs))
      res
    }
  )
)

#' Class for fitting a node sampling
nodeSampling_fit <-
  R6::R6Class(classname = "nodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param ... used for compatibility
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "node")
      private$psi <- private$card_N_o / (private$card_N_o + private$card_N_m)
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      res <- private$card_N_o * log(private$psi + 1 * (private$psi == 0)) +
        private$card_N_m * log(1 - private$psi + 1 * (private$psi == 1) )
      res
    }
  )
)

#' Class for fitting a node-centered sampling with covariate
covarNodeSampling_fit <-
  R6::R6Class(classname = "covarNodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  #' @description constructor
  #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
  #' @param ... used for compatibility
  public = list(
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "covar-node")
      y <- 1 * (partlyObservedNetwork$observedNodes)
      glm_out     <- suppressWarnings(glm.fit(cbind(1, partlyObservedNetwork$covarMatrix), y, family = binomial()))
      private$psi <- coefficients(glm_out)
      private$rho <- fitted(glm_out)
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      res <- sum(log(private$rho[private$N_obs])) + sum(log(1 - private$rho[!private$N_obs]))
      res
    }
  )
)

#' Class for fitting a double-standard sampling
#'
#' A dyad is observed with probability \eqn{p_0} if it is 0, \eqn{p_1} if it is 1 (the vector of
#' parameters is \eqn{(p_0, p_1)}). Given the blocks, a dyad between blocks q and l is missing
#' with probability \eqn{m_{ql} = \pi_{ql}(1 - p_1) + (1 - \pi_{ql})(1 - p_0)}: the missing dyads
#' are integrated out given the blocks, and the posterior probability that a missing dyad between
#' blocks q and l is an edge is \eqn{\nu_{ql} = \pi_{ql}(1 - p_1) / m_{ql}}, a Q x Q matrix.
#' The parameters are updated jointly with the connectivity of the SBM (see
#' [SimpleSBM_fit_MNAR]): this class only stores them, with the terms of the sampling.
#'
#' With a single block, \eqn{(\pi, p_0, p_1)} are not identifiable: the observed data only inform
#' \eqn{\pi p_1} and \eqn{(1 - \pi) p_0}. The fit is then constrained to \eqn{p_0 = p_1}, that is,
#' to dyad sampling, with one parameter.
doubleStandardSampling_fit <-
  R6::R6Class(classname = "doubleStandardSampling_fit",
  inherit = networkSamplingDyads_fit,
  private = list(
    So     = NULL, ## number of observed edges
    So.bar = NULL, ## number of observed non-edges
    S      = NULL, ## pattern of the missing dyads (upper triangle if undirected)
    Z      = NULL, ## memberships of the last update of the parameters
    ZtSZ   = NULL, ## Z' S Z: number of missing dyads between blocks q and l
    theta  = NULL, ## connectivity of the SBM at the last update of the parameters
    ## probability that a dyad between blocks q and l is missing
    m      = function() private$theta * (1 - private$psi[2]) + (1 - private$theta) * (1 - private$psi[1]),
    ## posterior probability that a missing dyad between blocks q and l is an edge
    nu     = function() check_boundaries(private$theta * (1 - private$psi[2]) / private$m()),
    single_block = function() ncol(private$Z) == 1,
    set_blocks = function(Z, theta) {
      private$Z     <- Z
      private$ZtSZ  <- as.matrix(t(Z) %*% private$S %*% Z)
      private$theta <- theta
    }
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    #' @param connectInit Q x Q matrix of initial connectivity of the SBM
    initialize = function(partlyObservedNetwork, blockInit, connectInit) {
      super$initialize(partlyObservedNetwork, "double-standard")
      private$So     <- sum(partlyObservedNetwork$networkData)
      private$So.bar <- private$card_D_o - private$So
      private$S      <- partlyObservedNetwork$samplingMatrixBar
      ## starts from dyad sampling, p0 = p1
      private$psi    <- rep(private$card_D_o / private$card_D, 2)
      private$set_blocks(blockInit, connectInit)
    },
    #' @description stores the parameters, updated jointly with the connectivity of the SBM
    #'   (see [SimpleSBM_fit_MNAR]), and the current blocks
    #' @param Z probabilities of block memberships
    #' @param theta the connectivity of the SBM
    #' @param psi the updated (p0, p1); \code{NULL} with a single block, where p0 = p1 is the
    #'   observation rate (dyad sampling)
    #' @param ... use for compatibility
    update_parameters = function(Z, theta, psi = NULL, ...) {
      private$set_blocks(Z, theta)
      private$psi <- if (private$single_block() || is.null(psi))
        rep(private$card_D_o / private$card_D, 2) else psi
      invisible(self)
    },
    #' @description no imputation is needed along the VEM: the missing dyads are integrated out
    #'   given the blocks (see the field \code{blockImputation})
    #' @param ... use for compatibility
    #' @return \code{NULL}
    update_imputation = function(...) {NULL},
    #' @description permute the parameters indexed by the blocks, after the blocks of the SBM
    #'   were permuted
    #' @param o the permutation applied to the blocks (new block k is former block \code{o[k]})
    reorder = function(o) {
      private$Z     <- private$Z[, o, drop = FALSE]
      private$ZtSZ  <- private$ZtSZ[o, o, drop = FALSE]
      private$theta <- private$theta[o, o, drop = FALSE]
      invisible(self)
    }
  ),
  active = list(
    #' @field df the number of parameters: 2, or 1 with a single block (p0 = p1)
    df = function(value) {if (private$single_block()) 1 else 2},
    #' @field vExpec variational expectation of the log-likelihood of the sampling, the missing
    #'   dyads being integrated out given the blocks:
    #'   \eqn{E_\tau[\log p(R \mid Y^o, Z)]}
    vExpec = function(value) {
      log(private$psi[2]) * private$So + log(private$psi[1]) * private$So.bar +
        sum(private$ZtSZ * log(private$m()))
    },
    #' @field log_lambda Q x Q matrix, the log-probability that a dyad between blocks q and l is
    #'   missing, the term of the missing dyads in the update of the blocks
    log_lambda = function(value) {log(private$m())},
    #' @field blockImputation Q x Q matrix, the posterior probability that a missing dyad between
    #'   blocks q and l is an edge, with the current parameters
    blockImputation = function(value) {private$nu()}
  )
)

#' Class for fitting a block-dyad sampling
blockDyadSampling_fit <-
  R6::R6Class(classname = "blockDyadSampling_fit",
  inherit  = networkSamplingDyads_fit,
  private  = list(
    R        = NULL,  ## sampling matrix
    S        = NULL,
    Z        = NULL,
    directed = NULL
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(partlyObservedNetwork, blockInit) {
      super$initialize(partlyObservedNetwork, "block-dyad")
      private$R <- partlyObservedNetwork$samplingMatrix
      private$S <- partlyObservedNetwork$samplingMatrixBar
      private$directed <- partlyObservedNetwork$is_directed
      self$update_parameters(NA, blockInit)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    #' @param Z probabilities of block memberships
    update_parameters = function(nu, Z) {
      private$Z <- Z
      ZtRZ <- as.matrix(t(Z) %*% private$R %*% Z)
      Zbar <- colSums(Z)
      if(private$directed) {
        private$psi <- check_boundaries(ZtRZ / ( Zbar %o% Zbar - Zbar ))
      } else {
        private$psi <- check_boundaries(( ZtRZ + t(ZtRZ) ) / ( Zbar %o% Zbar - Zbar ))
      }
    },
    #' @description permute the parameters indexed by the blocks, after the blocks of the SBM
    #'   were permuted
    #' @param o the permutation applied to the blocks (new block k is former block \code{o[k]})
    reorder = function(o) {
      private$psi <- private$psi[o, o, drop = FALSE]
      private$Z   <- private$Z[, o, drop = FALSE]
      invisible(self)
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {
      ## E_q[log p(R | Z)] = sum_{(i,j) observed} sum_{q,l} tau_iq tau_jl log(psi_ql) + same on the
      ## missing dyads with log(1 - psi_ql): the expectation of the log, consistent with the update
      ## of psi and with log_lambda (it was the log of the expectation, log(Z psi Z^T), before), and
      ## computed through Z^T R Z and Z^T S Z, without any n x n matrix
      Z <- private$Z
      sum(as.matrix(t(Z) %*% private$R %*% Z) * log(private$psi)) +
        sum(as.matrix(t(Z) %*% private$S %*% Z) * log(1 - private$psi))
    },
    #' @field log_lambda matrix, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {
      res <- private$R %*% private$Z %*% t(log(private$psi)) + private$S %*% private$Z %*% t(log(1 - private$psi)) +
        t(private$R) %*% private$Z %*% log(private$psi) + t(private$S) %*% private$Z %*% log(1 - private$psi)
      res
    }
  )
)

#' Class for fitting a block-node sampling
blockNodeSampling_fit <-
  R6::R6Class(classname = "blockNodeSampling_fit",
  inherit = networkSamplingNodes_fit,
  private = list(
    So     = NULL, ## sum_(i in Nobs ) Z_iq
    Sm     = NULL  ## sum_(i in Nmiss) Z_iq
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(partlyObservedNetwork, blockInit) {
      super$initialize(partlyObservedNetwork, "block-node")
      self$update_parameters(NA, blockInit)
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param imputedNet an adjacency matrix where missing values have been imputed
    #' @param Z indicator of blocks
    update_parameters = function(imputedNet, Z) {
      private$So  <- colSums(Z[ private$N_obs, , drop = FALSE])
      private$Sm  <- colSums(Z[!private$N_obs, , drop = FALSE])
      private$psi <- check_boundaries(private$So / (private$So + private$Sm))
    },
    #' @description permute the parameters indexed by the blocks, after the blocks of the SBM
    #'   were permuted
    #' @param o the permutation applied to the blocks (new block k is former block \code{o[k]})
    reorder = function(o) {
      private$psi <- private$psi[o]
      private$So  <- private$So[o]
      private$Sm  <- private$Sm[o]
      invisible(self)
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      res <- c(crossprod(private$So, log(private$psi)) +  crossprod(private$Sm, log(1 - private$psi)))
      res
    },
    #' @field log_lambda double, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {
      res <- sapply(private$psi, function(psi) ifelse(private$N_obs, log(psi), log(1 - psi)))
      res
    }
  )
)

#' Class for fitting a degree sampling
degreeSampling_fit <-
  R6::R6Class(classname = "degreeSampling_fit",
  inherit = networkSamplingNodes_fit,
  ## Everything is stored in O(n + number of missing dyads): the degrees as vectors, the
  ## imputations on the sparse pattern of the missing dyads (upper triangle if undirected).
  private = list(
    directed = NULL,
    D_obs = NULL, # degrees on the observed part of the network (constant)
    miss_i = NULL, miss_j = NULL, # row and column of each missing dyad, in the order of the pattern
    miss_t = NULL, # directed: position of dyad (j, i) among the missing dyads, NA if it is observed
    nu_x = NULL, # current imputations of the missing dyads, in the order of the pattern
    ksi = NULL, # additional variational parameter for approximation in the logistic function
    D   = NULL  # expected degrees: observed part + imputations of the missing dyads
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param blockInit n x Q matrix of initial block indicators
    #' @param connectInit Q x Q matrix of initial block probabilities of connection
    initialize = function(partlyObservedNetwork, blockInit, connectInit) {
      super$initialize(partlyObservedNetwork, "degree")

      private$directed <- partlyObservedNetwork$is_directed
      S <- partlyObservedNetwork$samplingMatrixBar
      private$miss_i <- S@i + 1L
      private$miss_j <- rep.int(seq_len(ncol(S)), diff(S@p))
      if (private$directed) {
        n <- ncol(S)
        key <- private$miss_i + (private$miss_j - 1) * n
        private$miss_t <- match(private$miss_j + (private$miss_i - 1) * n, key)
      }
      Y <- partlyObservedNetwork$networkData
      private$D_obs <- as.numeric(
        if (private$directed) Matrix::rowSums(Y) else Matrix::rowSums(Y) + Matrix::colSums(Y)
      )
      ## degrees on the observed part, to initialize psi
      private$D <- private$D_obs

      ## will fluctuate along the algorithm
      private$psi <- coefficients(glm(1*(private$N_obs) ~ private$D, family = binomial(link = "logit")))

      self$update_parameters(.bilinear_at_pattern(blockInit, connectInit, S))
    },
    #' @description a method to update the estimation of the parameters. By default, nothing to do (corresponds to MAR sampling)
    #' @param nu the imputations of the missing dyads (sparse matrix with the pattern of the missing dyads)
    #' @param ... used for compatibility
    ## Jaakkola-Jordan variational bound on the logit(psi[1] + psi[2]*degree) sampling model:
    ## ksi is the per-node bound parameter, h(ksi) its associated quadratic-bound coefficient
    ## (see utils_missSBM.R), (a,b) solve the resulting 2x2 weighted least-squares system for
    ## psi -- see Tabouy, Barbillon & Chiquet (2019), \doi{10.1080/01621459.2018.1562934}
    update_parameters = function(nu, ...) {
      ## nu: imputations of the missing dyads, on the pattern of the missing dyads
      private$nu_x <- nu@x
      ## expected degrees and their second moments: sums over the rows (and the columns, for the
      ## other triangle of an undirected network) of the imputations
      row_sums <- function(M) {
        as.numeric(if (private$directed) Matrix::rowSums(M) else Matrix::rowSums(M) + Matrix::colSums(M))
      }
      private$D <- private$D_obs + row_sums(nu)
      V <- nu
      V@x <- nu@x * (1 - nu@x)
      D2 <- row_sums(V) + private$D^2
      private$ksi <- check_boundaries(sqrt( private$psi[1]^2 + private$psi[2]^2 * D2  + 2 * private$psi[1] * private$psi[2] * private$D))
      C  <- .5 * length(private$D) - sum(!private$N_obs)
      s1 <- .5 * sum(private$D) - sum(private$D[!private$N_obs])
      s_hksi   <- sum(h(private$ksi))
      s_hksiD  <- sum(h(private$ksi) * private$D)
      s_hksiD2 <- sum(h(private$ksi) * D2)
      denom <- 2 * (s_hksi * s_hksiD2 - s_hksiD^2)
      a <- (s1 * s_hksiD - C * s_hksiD2) / denom
      b <- (C * s_hksiD - s1 * s_hksi) / denom
      private$psi <- c(a, b)
    },
    #' @description a method to update the imputation of the missing entries.
    #' @param PI the matrix of inter/intra class probability of connection
    #' @param ... use for compatibility
    update_imputation = function(PI,...) {
      ## for a missing dyad (i, j): C_ij + C_ji - psi[2], with
      ## C_ij = 2 h(ksi_i) (psi[1] psi[2] + psi[2]^2 (1 + D_i - nu_ij)), nu_ij taken as 0 when the
      ## dyad is observed (only possible for (j, i) in a directed network)
      i <- private$miss_i; j <- private$miss_j
      nu_ij <- private$nu_x
      nu_ji <- if (private$directed) {
        res <- private$nu_x[private$miss_t]; res[is.na(res)] <- 0; res
      } else nu_ij
      hk <- 2 * h(private$ksi)
      a  <- private$psi[1] * private$psi[2]
      b  <- private$psi[2]^2
      eta <- hk[i] * (a + b * (1 + private$D[i] - nu_ij)) + hk[j] * (a + b * (1 + private$D[j] - nu_ji)) - private$psi[2]
      nu  <- PI
      nu@x <- check_boundaries(.logistic(.logit(PI@x) + eta))
      nu
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function() {
      prob <-  check_boundaries(.logistic(private$psi[1] + private$psi[2] * private$D))
      res  <-  sum(private$N_obs * log(prob)) + sum( (!private$N_obs) * log(1 - prob))
      res
    }
  )
)
