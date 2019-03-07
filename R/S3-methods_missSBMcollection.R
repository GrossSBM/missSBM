#' @rdname missSBM
#' @export
is.missSBMcollection <- function(Robject) {
  inherits(Robject, "missSBMcollection")
}

#' @rdname missSBM
#' @export
ICL <- function(Robject) { UseMethod("ICL", Robject) }

#' @rdname missSBM
#' @importFrom stats setNames
#' @export
ICL.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  setNames(sapply(Robject, function(model) model$vICL), names(Robject))
}

#' @rdname missSBM
#' @export
optimizationStatus <- function(Robject) { UseMethod("optimizationStatus", Robject) }

#' @rdname missSBM
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

#' @rdname missSBM
#' @export
getBestModel <- function(Robject) {UseMethod("getBestModel", Robject)}

#' @rdname missSBM
#' @export
getBestModel.missSBMcollection <- function(Robject) {
  stopifnot(is.missSBMcollection(Robject))
  Robject[[which.min(ICL(Robject))]]
}

