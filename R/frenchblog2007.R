#' Political Blogosphere network prior to 2007 French presidential election
#'
#' French Political Blogosphere network dataset consists of a single day snapshot of
#' over 200 political blogs automatically extracted the 14 October 2006 and manually classified by
#' the "Observatoire Presidentielle" project. Originally part of the 'mixer' package
#'
#' @format An igraph object with 196 nodes. The vertex attribute "party" provides a possible
#' clustering of the nodes.
#' @source \url{https://www.linkfluence.com/}
#' @examples
#' data(frenchblog2007)
#' igraph::V(frenchblog2007)$party
#' igraph::plot.igraph(frenchblog2007,
#'   vertex.color = factor(igraph::V(frenchblog2007)$party),
#'   vertex.label = NA
#'  )
"frenchblog2007"
