#' Virtual class of the fits of a sampling design
#'
#' The mother of [`networkSamplingDyads_fit`] and [`networkSamplingNodes_fit`], with the defaults
#' of a MAR design: no parameter indexed by the blocks, nothing to update along the VEM.
#' @include R6Class-networkSampling.R
#' @import R6
networkSampling_fit <-
  R6::R6Class(classname = "networkSampling_fit",
  inherit = networkSampling,
  private = list(
    card = NULL # number of sampling units (dyads or nodes), counted by the penalty
  ),
  public = list(
    #' @description show method
    show = function() {
      super$show()
      cat("  $penalty, $vExpec\n")
    },
    #' @description update the parameters along the VEM: nothing to do by default, the
    #'   parameters of a MAR design not depending on the SBM
    #' @param ... use for compatibility
    update_parameters = function(...) {invisible(NULL)},
    #' @description a method to update the imputation of the missing entries.
    #' @param nu the matrix of (uncorrected) imputation for missing entries
    update_imputation = function(nu) {nu},
    #' @description permute the parameters indexed by the blocks, after the blocks of the SBM
    #'   were permuted. By default, nothing to do (no parameter indexed by the blocks)
    #' @param o the permutation applied to the blocks (new block k is former block \code{o[k]})
    reorder = function(o) {invisible(self)},
    #' @description the variational expectation of the log-likelihood of the sampling at given
    #'   memberships (used by the ICL, at the MAP memberships)
    #' @param Z an n x Q matrix of memberships
    vExpec_at = function(Z) {self$vExpec} # the sampling does not depend on the blocks
  ),
  active = list(
    #' @field penalty double, value of the penalty term in ICL: the parameters of the sampling
    #'   count all the sampling units (dyads or nodes), the sampling matrix being fully observed
    penalty = function(value) {log(private$card) * self$df},
    #' @field log_lambda the term of the sampling in the update of the blocks: 0 when it does not
    #'   depend on them
    log_lambda = function(value) {0}
  )
)

#' Virtual class of the fits of the sampling designs on the dyads
#' @import R6
networkSamplingDyads_fit <-
  R6::R6Class(classname = "networkSamplingDyads_fit",
  inherit = networkSampling_fit,
  private = list(
    card_D   = NULL, # number of possible dyads in the network
    card_D_o = NULL, # number of observed dyads in the network
    card_D_m = NULL  # number of missing dyads in the network
  ),
  public = list(
    #' @description constructor
    #' @param partlyObservedNetwork a object with class partlyObservedNetwork representing the observed data with possibly missing entries
    #' @param name a character for the name of sampling to fit on the partlyObservedNetwork
    initialize = function(partlyObservedNetwork, name) {
      private$name     <- name
      private$card_D   <- partlyObservedNetwork$nbDyads
      private$card_D_o <- sum(partlyObservedNetwork$samplingMatrix)
      private$card_D_m <- private$card_D - private$card_D_o
      private$card     <- private$card_D
    }
  )
)

#' Virtual class of the fits of the sampling designs on the nodes
#' @import R6
networkSamplingNodes_fit <-
  R6::R6Class(classname = "networkSamplingNodes_fit",
  inherit = networkSampling_fit,
  private = list(
    card_N   = NULL, # number of nodes in the network
    card_N_o = NULL, # number of observed nodes
    card_N_m = NULL, # number of missing nodes
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
      private$card     <- private$card_N
    }
  )
)

#' Class for fitting a dyad sampling
dyadSampling_fit <-
  R6::R6Class(classname = "dyadSampling_fit",
  inherit = networkSamplingDyads_fit,
  public = list(
    #' @description constructor
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
    #' @param ... used for compatibility
    initialize = function(partlyObservedNetwork, ...) {
      super$initialize(partlyObservedNetwork, "covar-dyad")
      dyads <- which(valid_dyads(dim(partlyObservedNetwork$networkData), partlyObservedNetwork$is_directed), arr.ind = TRUE)
      X <- cbind(1, apply(partlyObservedNetwork$covarArray, 3, function(x) x[dyads]))
      y <- partlyObservedNetwork$samplingMatrix[dyads]
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
  public = list(
    #' @description constructor
    #' @param ... used for compatibility
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
    ## the log-likelihood of the sampling given the numbers of missing dyads between the blocks
    vExpec_given = function(ZtSZ) {
      log(private$psi[2]) * private$So + log(private$psi[1]) * private$So.bar + sum(ZtSZ * log(private$m()))
    },
    set_blocks = function(Z, theta) {
      private$Z     <- Z
      private$ZtSZ  <- .block_counts(Z, private$S)
      private$theta <- theta
    }
  ),
  public = list(
    #' @description constructor
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
    },
    #' @description the variational expectation of the log-likelihood of the sampling at given
    #'   memberships (used by the ICL, at the MAP memberships)
    #' @param Z an n x Q matrix of memberships
    vExpec_at = function(Z) {private$vExpec_given(.block_counts(Z, private$S))}
  ),
  active = list(
    #' @field df the number of parameters: 2, or 1 with a single block (p0 = p1)
    df = function(value) {if (private$single_block()) 1 else 2},
    #' @field vExpec variational expectation of the log-likelihood of the sampling, the missing
    #'   dyads being integrated out given the blocks:
    #'   \eqn{E_\tau[\log p(R \mid Y^o, Z)]}
    vExpec = function(value) {private$vExpec_given(private$ZtSZ)},
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
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(partlyObservedNetwork, blockInit) {
      super$initialize(partlyObservedNetwork, "block-dyad")
      private$R <- partlyObservedNetwork$samplingMatrix
      private$S <- partlyObservedNetwork$samplingMatrixBar
      private$directed <- partlyObservedNetwork$is_directed
      self$update_parameters(blockInit)
    },
    #' @description update the parameters: the observation rate of the dyads between each pair
    #'   of blocks
    #' @param Z probabilities of block memberships
    #' @param ... use for compatibility
    update_parameters = function(Z, ...) {
      private$Z <- Z
      ZtRZ <- .block_counts(Z, private$R)
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
    },
    #' @description the variational expectation of the log-likelihood of the sampling at given
    #'   memberships (used by the ICL, at the MAP memberships)
    #' @param Z an n x Q matrix of memberships
    vExpec_at = function(Z) {
      ## E_q[log p(R | Z)] = sum_{(i,j) observed} sum_{q,l} tau_iq tau_jl log(psi_ql) + same on the
      ## missing dyads with log(1 - psi_ql): the expectation of the log (not the log of the
      ## expectation), consistent with the update of psi and with log_lambda, computed through
      ## Z^T R Z and Z^T S Z, without any n x n matrix
      sum(.block_counts(Z, private$R) * log(private$psi)) +
        sum(.block_counts(Z, private$S) * log(1 - private$psi))
    }
  ),
  active = list(
    #' @field vExpec variational expectation of the sampling
    vExpec = function(value) {self$vExpec_at(private$Z)},
    #' @field log_lambda matrix, term for adjusting the imputation step which depends on the type of sampling
    log_lambda = function(value) {
      ## crossprod(R, Z) = t(R) %*% Z without forming the transpose of R
      res <- private$R %*% private$Z %*% t(log(private$psi)) + private$S %*% private$Z %*% t(log(1 - private$psi)) +
        Matrix::crossprod(private$R, private$Z) %*% log(private$psi) + Matrix::crossprod(private$S, private$Z) %*% log(1 - private$psi)
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
    #' @param blockInit n x Q matrix of initial block indicators
    initialize = function(partlyObservedNetwork, blockInit) {
      super$initialize(partlyObservedNetwork, "block-node")
      self$update_parameters(blockInit)
    },
    #' @description update the parameters: the observation rate of the nodes of each block
    #' @param Z probabilities of block memberships
    #' @param ... use for compatibility
    update_parameters = function(Z, ...) {
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
    },
    #' @description the variational expectation of the log-likelihood of the sampling at given
    #'   memberships (used by the ICL, at the MAP memberships)
    #' @param Z an n x Q matrix of memberships
    vExpec_at = function(Z) {
      c(crossprod(colSums(Z[ private$N_obs, , drop = FALSE]), log(private$psi)) +
        crossprod(colSums(Z[!private$N_obs, , drop = FALSE]), log(1 - private$psi)))
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
    D   = NULL, # expected degrees: observed part + imputations of the missing dyads
    D2  = NULL  # expected squared degrees
  ),
  public = list(
    #' @description constructor
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
    #' @description update the parameters given the imputations: the expected degrees, the
    #'   parameters zeta of the bound of Jaakkola and Jordan, then (a, b)
    #' @param nu the imputations of the missing dyads (sparse matrix with the pattern of the missing dyads)
    #' @param ... used for compatibility
    update_parameters = function(nu, ...) {
      ## bound of Jaakkola and Jordan on the logistic model logit(psi[1] + psi[2] * degree): zeta
      ## is the parameter of the bound of each node, h(zeta) its quadratic coefficient, and (a, b)
      ## solve the resulting 2 x 2 weighted least-squares system (Tabouy, Barbillon & Chiquet, 2019)
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
      private$D2 <- D2
      ## the optimal zeta_i = sqrt(E[(a + b D_i)^2]), any positive value (not a probability: it is
      ## often above 1); bounded away from 0, where h(zeta) is 0/0
      private$ksi <- pmax(sqrt(pmax(private$psi[1]^2 + private$psi[2]^2 * D2 + 2 * private$psi[1] * private$psi[2] * private$D, 0)), 1e-6)
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
      ## in one pass over the missing dyads, in C++
      nu  <- PI
      nu@x <- degree_imputation_cpp(PI@x, private$miss_i, private$miss_j, private$nu_x,
                                    if (private$directed) private$miss_t else integer(0),
                                    2 * h(private$ksi), private$D, private$psi[1], private$psi[2],
                                    private$directed)
      nu
    }
  ),
  active = list(
    #' @field vExpec lower bound of the variational expectation of the log-likelihood of the
    #'   sampling, by the bound of Jaakkola and Jordan: with x_i = a + b D_i,
    #'   \eqn{\log \sigma(x) \geq \log \sigma(\zeta) + (x - \zeta)/2 + h(\zeta)(x^2 - \zeta^2)}
    #'   and \eqn{\log(1 - \sigma(x)) = \log \sigma(x) - x}, in expectation over the imputations
    vExpec = function() {
      a <- private$psi[1]; b <- private$psi[2]; z <- private$ksi; D <- private$D
      - sum((a + b * D)[!private$N_obs]) +
        sum(-log1p(exp(-z)) + (a + b * D - z) / 2 + h(z) * (a^2 + 2 * a * b * D + b^2 * private$D2 - z^2))
    }
  )
)
