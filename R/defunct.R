#' @name missSBM-defunct
#' @title Defunct Functions in Package missSBM
NULL

#' @export
#' @rdname missSBM-defunct
#' @details `estimate()` is replaced by `missSBM::estimateMissSBM`.
#' @param ... unused arguments
estimate <- function(...) {
  .Defunct("missSBM::estimateMissSBM")
}

#' @export
#' @rdname missSBM-defunct
#' @details `sample()` is replaced by `missSBM::observeNetwork`.
#' @param ... unused arguments
sample <- function(...) {
  .Defunct("missSBM::observeNetwork")
}

#' @export
#' @rdname missSBM-defunct
#' @details `simualte()` is replaced by `sbm::sampleSimpleSBM`.
#' @param ... unused arguments
simulate <- function(...) {
  .Defunct("sbm::sampleSimpleSBM")
}

