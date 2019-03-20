#' An R package for adjusting Stochastic Block Models from networks data sampled under various missing data conditions
#'
#' The SBM package provides three main functions: \code{\link{simulate}}, \code{\link{sample}}, \code{\link{estimate}}.
#'
#' @docType package
#' @author Timothée Tabouy, Pierre Barbillon, Julien Chiquet
#' @references Timothée Tabouy, Pierre Barbillon & Julien Chiquet (2019) “Variational Inference for Stochastic Block Models from Sampled Data”, Journal of the American Statistical Association, <doi:10.1080/01621459.2018.1562934>
#' @import R6 methods
#' @useDynLib missSBM
#' @importFrom Rcpp sourceCpp
#' @name missSBM
NULL
