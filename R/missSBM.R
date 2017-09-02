#' @title fit a series of SBM on network data possible sampled under various missing data condition
#'
#' @field network an R object with class sampledNetwork describing the network data. See the corresponding documentation and examples
#' @field vBlocks a numeric vector containing the number of blocks considered for the successively adjusted SBM
#' @field family a character indicating the distribution family of the SBM
#' @field sampling a character indicating the modelign of the missing data in the sampled network. Default is MAR.
#' @field control a list to control the algorithm. See details.
#'
#' @export
missSBM <- function(sampledNetwork,
                    vBlocks,
                    family   = "Bernoulli",
                    sampling = "MAR",
                    control=list()) {
  ## TODO: everything
  ## will used all classes
}

possible.sampling <- c("MAR.edge", "MAR.node", "double.standard", "star.degree", "class")
possible.family   <- c("BernoulliDirected", "BernoulliUndirected", "PoissonDirected", "PoissonUndirected")

    # ## Check sampling design and adequency to the sampled data
    # if (sum(is.na(network)) == 0) {
    #   self$sampling <- "full"
    # } else {
    #   if (is.null(sampling)) {
    #     self$sampling <- "MAR.edge"
    #   } else {
    #     stopifnot(sampling %in% possible.sampling)
    #     self$sampling <- sampling
    #   }
    # }
    #
    # ## Check adequency of distribution choice and sampled data
    # stopifnot(all(network[!is.na(network)] >= 0))
    # if (range(network, na.rm=TRUE) == c(0,1)) {
    #   if (is.null(distribution)) {
    #     distribution <- "Bernoulli"
    #   } else {
    #     stopifnot(distribution == "Bernoulli")
    #   }
    # }
