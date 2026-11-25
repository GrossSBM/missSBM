# Summary method for a [`missSBM_collection`](https://grosssbm.github.io/missSBM/reference/missSBM_collection.md)

Summary method for a
[`missSBM_collection`](https://grosssbm.github.io/missSBM/reference/missSBM_collection.md)

## Usage

``` r
# S3 method for class 'missSBM_collection'
summary(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`missSBM_collection`](https://grosssbm.github.io/missSBM/reference/missSBM_collection.md)

- ...:

  additional parameters for S3 compatibility.

## Value

a basic printing output

## Examples

``` r
## a fit of the French political blogosphere, 30% of its dyads missing
adjMatrix <- missSBM::frenchblog2007 %>%
  igraph::as_adjacency_matrix(sparse = FALSE) %>%
  missSBM::observeNetwork(sampling = "dyad", parameters = 0.7)
collection <- estimateMissSBM(adjMatrix, 1:3, "dyad",
  control = missSBM_param(trace = FALSE, iterates = 0))
summary(collection)
#> --------------------------------------------------------
#> COLLECTION OF 3 SBM fits          
#> ========================================================
#>  - Number of blocks considered: from 1 to 3
#>  - Best model (smaller ICL): 3
#>  - Fields: $models, $ICL, $vBlocks, $occupiedBlocks, $degenerate, $bestModel, $optimizationStatus
#>  - Method: $estimate(), $polish(), $explore(), $plot() 
```
