zero <- .Machine$double.eps

available_samplings <- c("dyad", "node", "snowball", "degree", "block", "double_standard")

available_samplings_covariates <- c("dyad", "node")

bar <- function(X) {
  X.bar <- 1 - X ; diag(X.bar) <- 0
  X.bar
}

quad_form <- function(A,x) {t(x) %*% A %*% x}

logistic <- function(x) {1/(1 + exp(-x))}
logit    <- function(x) {log(x/(1 - x))}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

h <- function(x) {-.5 * (logistic(x) - 0.5) / x}

xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x*log(x))

check_boundaries <- function(x) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  x
}

#' @importFrom graphics box image par
#' @export
image_NA <- function(z,  zlim = c(0,1), col = c("white", "midnightblue"), na.color = 'gray', outside.below.color = 'black', outside.above.color = 'white', ...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA

  z[which(z < zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z > zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z > zlim[2]))] <- newz.na # same for newz.na

  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na

  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range

  par(mar = c(2.1,8.1,3.1,3.1))
  image(z[nrow(z):1,],  zlim = zlim, col = col, xaxt = "n", yaxt = "n", ...) # we finally call image(...)
  box()
  # par(mar=c(5.1,4.1,4.1,2.1))
}


#' #' @import ggplot2 igraph viridis
#' #' @importFrom dplyr inner_join mutate select arrange
#' #' @export
#' gg_image_NA <- function(adjacencyMatrix, memberships) {
#'
#'   adjacencyMatrix[is.na(adjacencyMatrix)] <- -1
#'   G <- graph_from_adjacency_matrix(adjacencyMatrix, weighted = TRUE)
#'   V(G)$membership <- memberships
#'   V(G)$name <- 1:ncol(adjacencyMatrix)
#'   E(G)$miss <- E(G)$weight == -1
#'   node_list <- get.data.frame(G, what = "vertices")
#'
#'   edge_list <- get.data.frame(G, what = "edges") %>%
#'     inner_join(node_list %>% select(name, membership), by = c("from" = "name")) %>%
#'     inner_join(node_list %>% select(name, membership), by = c("to" = "name")) %>%
#'     mutate(membership = ifelse(membership.x == membership.y, membership.x, "dyad") %>% factor()) %>%
#'     mutate(missing = ifelse(miss, "missing", "observed") %>% factor()) %>% select(-weight, -miss) %>%
#'     mutate(membership_missingness = paste(membership,missing, sep = "-"))
#'
#'   all_nodes <- sort(node_list$name)
#'
#'   plot_data <- edge_list %>% mutate(
#'     to = factor(to, levels = all_nodes),
#'     from = factor(from, levels = all_nodes))
#'   name_order <- (node_list %>% arrange(membership))$name
#'
#'   plot_data <- edge_list %>% mutate(
#'     to = factor(to, levels = name_order),
#'     from = factor(from, levels = name_order))
#'
#'   if (sum(plot_data$missing == "missing") == 0) {
#'     p <- ggplot(plot_data, aes(x = from, y = to, fill = membership))
#'   } else {
#'     p <- ggplot(plot_data, aes(x = from, y = to, fill = membership_missingness))
#'   }
#'   p <- p + geom_raster() + theme_classic() +
#'     scale_fill_viridis(discrete = TRUE, option = "magma") +
#'     theme(aspect.ratio = 1,
#'           axis.title = element_blank(),
#'           axis.text = element_blank() ,
#'           axis.ticks = element_blank(),
#'           axis.line = element_blank())
#'   p
#' }
