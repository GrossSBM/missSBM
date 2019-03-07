#' An R package for adjusting Stochastic Block Models from networks data sampled under various missing data conditions
#'
#' The SBM package provides three categories of important functions:
#' \code{\link{simulateSBM}}, \code{\link{sampleNetwork}}, \code{\link{missSBM}}.
#'
#' @docType package
#' @author Timothée Tabouy, Pierre Barbillon, Julien Chiquet
#' @references [1] Tabouy et al., Variationnal inference of Stochastic Block Model from sampled data (2017). arXiv:1707.04141.
#' @references [2] Kolaczyk Erik D. (2009). Statistical Analysis of Network Data: Methods and Models. Springer.
#' @import R6 methods
#' @useDynLib missSBM
#' @importFrom Rcpp sourceCpp
#' @name missSBM
NULL
