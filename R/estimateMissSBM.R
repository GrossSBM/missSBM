#' Control of a missSBM fit
#'
#' Helper to define the list of parameters that controls [estimateMissSBM()]. All arguments have defaults.
#'
#' @param threshold the VEM stops when an optimization step changes the objective function
#'   or the parameters by less than threshold. Default is 1e-2.
#' @param maxIter the VEM stops when the number of iterations exceeds maxIter. Default is 50.
#' @param fixPointIter number of fixed-point iterations in the variational E step. Default is 3.
#' @param imputation character, the imputation strategy used to build the initial clustering (see
#'   \code{\link{partlyObservedNetwork}}'s \code{imputation()}). Either "median" (default),
#'   "average" or "zero".
#' @param similarity an R x R -> R function to compute similarities between node covariates.
#'   Default is \code{\link{l1_similarity}}, that is, -abs(x - y). Only relevant when the
#'   covariates are node-centered (i.e. \code{covariates} is a list of size-N vectors).
#' @param useCov logical. If \code{covariates} is not empty, should they be used for the SBM
#'   inference (or just for the sampling)? Default is TRUE.
#' @param clusterInit initial clustering: \code{NULL} (default) for a cold spectral clustering, or
#'   a list with \code{length(vBlocks)} vectors, each with size \code{ncol(adjacencyMatrix)},
#'   providing a user-defined clustering.
#' @param polish logical, should each model be node-swap-polished (see [missSBM_fit]'s
#'   \code{polish()}) after the initial VEM fit, fixing individually misclassified nodes at that
#'   model's own number of blocks? Cheap relative to exploration (\code{iterates}), since it does
#'   not search across numbers of blocks. Default is TRUE.
#' @param iterates integer, the number of forward/backward exploration passes searching for a
#'   better number of blocks by splitting/merging blocks (unlike \code{polish}, which only
#'   refines each model at its own number of blocks): more expensive than \code{polish}.
#'   \code{0} disables exploration entirely (only the initial VEM fit, plus \code{polish} if
#'   enabled, is returned). Default is 1, enough on small networks; 2 are advised for a reference
#'   analysis of a network with many blocks.
#' @param maxMergeCandidates integer, caps the number of block-pair merge candidates tried
#'   during backward exploration (quadratic in the number of blocks otherwise). Beyond this cap,
#'   only the pairs with the most similar fitted connectivity profiles are tried, since merging
#'   two blocks with very different connectivity is rarely competitive anyway. Default is 30.
#' @param stopOnDegenerate logical. A requested number of blocks can be higher than what the
#'   network actually supports: VEM then collapses one or more classes (see [missSBM_fit]'s
#'   \code{repair()}), and even a fair recovery attempt (\code{repair()}, then a full exploration
#'   pass letting \code{explore_forward()} try to split a healthy neighbor into that range) can
#'   still collapse right back down. When this persists over \code{maxConsecutiveDegenerate}
#'   consecutive (increasing) values of \code{vBlocks} after a pass, forward (split) exploration
#'   stops growing further into that range on subsequent passes (only relevant when
#'   \code{iterates > 1}); the affected models stay as fitted, flagged via \code{$degenerate},
#'   and a warning is issued. Default is TRUE.
#' @param maxConsecutiveDegenerate integer, the run length that triggers \code{stopOnDegenerate}.
#'   Default is 2.
#' @param warmChain logical (\strong{work in progress}). If \code{TRUE}, each model is initialized
#'   by splitting (see [missSBM_fit]'s \code{split()}) the already-converged, smaller neighbor in
#'   \code{vBlocks} instead of an independent cold spectral clustering (see
#'   [missSBM_collection]'s \code{estimate_chain()}) -- meant to reduce VEM component collapse at
#'   higher numbers of blocks. Sequential in the number of blocks, unlike the default
#'   (\code{FALSE}), which fits every model in parallel: can be slower in wall-clock time when
#'   several workers are available. Default is FALSE.
#' @param trace logical for verbosity. Default is TRUE.
#'
#' @return a list of parameters configuring the fit, with class \code{missSBM_param}.
#'
#' @seealso [estimateMissSBM()]
#' @examples
#' my_control <- missSBM_param(iterates = 2, polish = FALSE)
#' my_control$iterates
#' @export
missSBM_param <- function(
    threshold                = 1e-2,
    maxIter                  = 50,
    fixPointIter             = 3,
    imputation                = c("median", "average", "zero"),
    similarity                = l1_similarity,
    useCov                    = TRUE,
    clusterInit               = NULL,
    polish                    = TRUE,
    iterates                  = 1,
    maxMergeCandidates        = 30,
    stopOnDegenerate          = TRUE,
    maxConsecutiveDegenerate  = 2,
    warmChain                 = FALSE,
    trace                     = TRUE
) {
  imputation <- match.arg(imputation)

  structure(list(
    threshold                = threshold,
    maxIter                  = maxIter,
    fixPointIter             = fixPointIter,
    imputation                = imputation,
    similarity                = similarity,
    useCov                    = useCov,
    clusterInit               = clusterInit,
    polish                    = polish,
    iterates                  = iterates,
    maxMergeCandidates        = maxMergeCandidates,
    stopOnDegenerate          = stopOnDegenerate,
    maxConsecutiveDegenerate  = maxConsecutiveDegenerate,
    warmChain                 = warmChain,
    trace                     = trace
  ), class = "missSBM_param")
}


#' Estimation of simple SBMs with missing data
#'
#' Variational EM inference of Stochastic Block Models indexed by block number from a partially observed network.
#'
#' @param adjacencyMatrix The N x N adjacency matrix of the network data. If \code{adjacencyMatrix} is symmetric,
#' we assume an undirected network with no loop; otherwise the network is assumed to be directed.
#' @param vBlocks The vector of number of blocks considered in the collection.
#' @param sampling The sampling design, the model of the process that originates the missing data:
#' MAR designs ("dyad", "node", "covar-dyad", "covar-node", "snowball") and MNAR designs
#' ("double-standard", "block-dyad", "block-node", "degree") are available. See details.
#' @param covariates An optional list with M entries (the M covariates). If the covariates are node-centered, each entry of \code{covariates}
#' must be a size-N vector;  if the covariates are dyad-centered, each entry of \code{covariates} must be N x N matrix.
#' @param control a list-like structure for controlling advanced features, with default generated
#' by \code{\link{missSBM_param}}. See the associated documentation for details.
#'
#' @return Returns an R6 object with class \code{\link{missSBM_collection}}.
#'
#' @details Internal functions use \code{future_lapply}, so set your plan to \code{'multisession'} or
#' \code{'multicore'} to use several cores/workers.
#'
#' @details The sampling designs are either Missing At Random (MAR) or Missing Not At Random
#' (MNAR), and observe either dyads or nodes (all the dyads of an observed node). See
#' \doi{10.1080/01621459.2018.1562934} for a complete description. The parameters of each design:
#' * Missing At Random (MAR)
#'    * "dyad": p, the probability that a dyad is observed.
#'    * "node": p, the probability that a node is observed.
#'    * "covar-dyad": \eqn{\beta \in R^M}, a dyad (i, j) being observed with probability
#'      \eqn{logistic(\beta' x_{ij})}, \eqn{x_{ij}} the covariates of the dyad.
#'    * "covar-node": \eqn{\beta \in R^M}, a node i being observed with probability
#'      \eqn{logistic(\beta' x_i)}, \eqn{x_i} the covariates of the node.
#'    * "snowball": (number of waves, proportion of the nodes observed in the first wave); the
#'      nodes of each wave are the neighbors of the nodes of the previous one. It is fitted as
#'      "node" sampling.
#' * Missing Not At Random (MNAR)
#'    * "double-standard": (p0, p1), the probabilities that a dyad is observed if it is 0 and if it
#'      is 1. (p0, p1) are identified through the differences of connectivity between the blocks:
#'      with a single block, they are not identifiable, and the model is fitted with p0 = p1, i.e.
#'      as dyad sampling (with a warning).
#'    * "block-node": \eqn{(p_1, \dots, p_Q)}, \eqn{p_q} the probability that a node of block q is
#'      observed.
#'    * "block-dyad": the Q x Q matrix \eqn{(p_{ql})}, \eqn{p_{ql}} the probability that a dyad
#'      between blocks q and l is observed.
#'    * "degree": (a, b), a node i being observed with probability \eqn{logistic(a + b D_i)},
#'      \eqn{D_i} its degree.
#'
#' @details The number of blocks is selected by the ICL (Integrated Classification Likelihood,
#' the smaller the better): -2 times the variational lower bound of the log-likelihood of the
#' observed data (the observed dyads and the sampling matrix), at the MAP memberships, plus a
#' penalty, in which the connectivity counts the observed dyads and the parameters of the sampling
#' all the dyads (or all the nodes). The ICLs of different designs fitted to the same data can be
#' compared to choose the design. When their samplings carry the same information (e.g.
#' "block-node" and "degree" on an assortative network, where the degree depends on the block),
#' the ICLs differ by a few units, about the penalty of one parameter (log N), and the partitions
#' are the same: such a difference is no evidence for one design. Each number of blocks is fitted
#' from a spectral clustering, then improved by a polish and by passes of exploration splitting
#' and merging blocks (see \code{\link{missSBM_param}}); one pass is enough on small networks,
#' two are advised for a reference analysis of a network with many blocks.
#'
#' @seealso \code{\link{observeNetwork}}, \code{\link{missSBM_collection}}, \code{\link{missSBM_fit}} and \code{\link{missSBM_param}}.
#' @examples
#' ## SBM parameters
#' N <- 100 # number of nodes
#' Q <- 3   # number of blocks
#' pi <- rep(1,Q)/Q     # block proportion
#' theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix
#'
#' ## Sampling parameters
#' samplingParameters <- .75 # the sampling rate
#' sampling  <- "dyad"      # the sampling design
#'
#' ## generate an undirected binary SBM with no covariate
#' sbm <- sbm::sampleSimpleSBM(N, pi, theta)
#'
#' ## Uncomment to set parallel computing with future
#' ## future::plan("multicore", workers = 2)
#'
#' ## Sample some dyads data + Infer SBM with missing data
#' collection <-
#'    observeNetwork(sbm$networkData, sampling, samplingParameters) %>%
#'    estimateMissSBM(vBlocks = 1:4, sampling = sampling)
#' plot(collection, "monitoring")
#' plot(collection, "icl")
#'
#' collection$ICL
#' coef(collection$bestModel$fittedSBM, "connectivity")
#'
#' myModel <- collection$bestModel
#' plot(myModel, "expected")
#' plot(myModel, "imputed")
#' plot(myModel, "meso")
#' coef(myModel, "sampling")
#' coef(myModel, "connectivity")
#' predict(myModel)[1:5, 1:5]
#'
#' @export
estimateMissSBM <- function(adjacencyMatrix, vBlocks, sampling, covariates = list(), control = missSBM_param()) {

  ## Temporary check for deprecated use of list()
  if (!inherits(control, "missSBM_param"))
    stop("We now use the function missSBM_param() to generate the list of parameters ",
         "that controls the fit: replace 'list(my_arg = xx)' by missSBM_param(my_arg = xx) ",
         "and see the documentation of missSBM_param().")

  ## Sanity checks
  stopifnot(sampling %in% available_samplings)
  stopifnot(is.numeric(vBlocks))
  stopifnot(is.character(sampling))
  stopifnot(is.list(covariates))

  ## explore_forward()/explore_backward()/estimate_chain() all assume vBlocks is strictly
  ## increasing (indices k, k +/- 1 refer to smaller/larger numbers of blocks)
  vBlocks_sorted <- sort(unique(vBlocks))
  if (length(vBlocks_sorted) != length(vBlocks) || any(vBlocks_sorted != vBlocks)) {
    warning("vBlocks was not strictly increasing: sorted and de-duplicated.", call. = FALSE)
  }
  vBlocks <- vBlocks_sorted

  ## If no covariate is provided, you cannot ask for using them
  if (length(covariates) == 0) control$useCov <- FALSE
  if (control$useCov) stopifnot(sampling %in% available_samplings_covariates)

  if (sampling == "double-standard" && 1 %in% vBlocks)
    warning("With a single block, double-standard sampling is not identifiable: ",
            "it is fitted with p0 = p1, i.e. as dyad sampling.", call. = FALSE)

  ## Prepare network data for estimation with missing data
  partlyObservedNet <- partlyObservedNetwork$new(adjacencyMatrix, covariates, control$similarity)
  clusterInit <- control$clusterInit
  if (is.null(clusterInit)) clusterInit <- partlyObservedNet$clustering(vBlocks, control$imputation)

  ## Instantiate the collection of missSBM_fit
  myCollection <- missSBM_collection$new(
      partlyObservedNet  = partlyObservedNet,
      sampling           = sampling,
      clusterInit        = clusterInit,
      control            = control
  )

  ## Launch estimation of each missSBM_fit
  if (isTRUE(control$warmChain)) myCollection$estimate_chain() else myCollection$estimate()

  ## Fix individually misclassified nodes at each model's own number of blocks
  if (control$polish) myCollection$polish()

  ## Looking for better models around (across numbers of blocks)
  myCollection$explore()

  ## Return the collection of adjusted missSBM_fit
  myCollection
}
