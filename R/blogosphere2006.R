#' French Political Blogosphere network
#'
#' French Political Blogosphere network dataset consists of a single day snapshot of
#' over 200 political blogs automatically extracted the 14 october 2006 and manually classified by
#' the ”Observatoire Presidentielle” project. Originally part of the 'mixer' package
#'
#' @format An igraph object with 196,. The vertex attribute "party" provides a possible clustering of the nodes.
#' @source \url{https://www.linkfluence.com/}
#' @examples
#' data(blogosphere2006)
#' V(blogosphere2006)
#' \dontrun{
#' igraph::plot(blogosphere2006)
#' }
#' @export
"blogosphere2006"
