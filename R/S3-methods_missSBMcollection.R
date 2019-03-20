#' @rdname estimate
#' @export
is.missSBMcollection <- function(Robject) {
  inherits(Robject, "missSBMcollection")
}

#' @rdname estimate
#' @export
ICL <- function(Robject) { UseMethod("ICL", Robject) }

#' @rdname estimate
#' @importFrom stats setNames
#' @export
ICL.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  setNames(sapply(Robject, function(model) model$vICL), names(Robject))
}

#' @rdname estimate
#' @export
optimizationStatus <- function(Robject) { UseMethod("optimizationStatus", Robject) }

#' @rdname estimate
#' @export
optimizationStatus.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  Reduce("rbind",
    lapply(Robject,
      function(model) {
        res <- model$monitoring
        res$nBlock <- model$fittedSBM$nBlocks
        res
    })
  )
}

#' @rdname estimate
#' @export
getBestModel <- function(Robject) {UseMethod("getBestModel", Robject)}

#' @rdname estimate
#' @export
getBestModel.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  Robject[[which.min(ICL(Robject))]]
}

