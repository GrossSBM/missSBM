library(tidyverse)
library(Matrix)
library(igraph)
library(rsvd)
library(missSBM)
library(ClusterR)

## =========================================================================
##
## Downloading data if required, formating into an igraph object

if(!file.exists("user_rating.txt.gz"))
  download.file("http://www.trustlet.org/datasets/extended_epinions/user_rating.txt.gz",
                "user_rating.txt.gz")
epinions_graph <- read_delim(
  "user_rating.txt.gz",
  delim = "\t",
  col_names = c("V1", "V2", "weight", "date")
  ) %>% dplyr::select(-4) %>% igraph::graph_from_data_frame(directed = TRUE)

## plot(degree_distribution(epinions), log = "xy")
min_degree <- 1
epinions_graph <- delete.vertices(epinions_graph, which(degree(epinions_graph) < min_degree))

## =========================================================================
##
## Initial clustering: 

L <- graph.laplacian(
  as.undirected(epinions_graph),
  normalized = TRUE,
  weights = NULL
 )
L@x[is.na(L@x) | is.infinite(L@x)] <- 0
svdL <- rsvd(L, nv = 20, nu = 20, k = 100)

v_nBlocks <- c(2, 3, 4, 5, 10, 15, 20)

clustering_init <- map(v_nBlocks, function(nBlocks) {
 km_model <-
   ClusterR::MiniBatchKmeans(
     data     = svdL$u[, 1:nBlocks, drop = FALSE],
     clusters = nBlocks, batch_size = 1000, num_init = 100
  )
 cl <- predict_MBatchKMeans(svdL$u[, 1:nBlocks, drop = FALSE], km_model$centroids)
 cl
}) %>% setNames(paste0(v_nBlocks, " blocks"))

counts <- map(clustering_init, tabulate)

A <- as_adjacency_matrix(epinions_graph, attr= "weight")
observed  <- as(1*(A !=  0), "dgCMatrix")
trusted   <- as(1*(A ==  1), "dgCMatrix")
untrusted <- as(1*(A == -1), "dgCMatrix")
epinions <- list(
  observed = as(1*(A !=  0), "dgCMatrix"),
  trusted  = as(1*(A ==  1), "dgCMatrix")
  )
