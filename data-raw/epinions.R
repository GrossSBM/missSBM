library(tidyverse)
library(igraph)

if(!file.exists("user_rating.txt.gz"))
  download.file("http://www.trustlet.org/datasets/extended_epinions/user_rating.txt.gz",
                "user_rating.txt.gz")

epinions_graph <- read_delim(
  "user_rating.txt.gz",
  delim = "\t",
  col_names = c("V1", "V2", "weight", "date")
  ) %>% select(-4) %>% igraph::graph_from_data_frame(directed = TRUE)

## plot(degree_distribution(epinions), log = "xy")
min_degree <- 1
epinions_subgraph <- delete.vertices(epinions_graph, which(degree(epinions_graph) < min_degree))

usethis::use_data(epinions_subgraph, overwrite = TRUE)
