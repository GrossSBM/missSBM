# Extract model fitted values from object [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md), return by [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)

Extract model fitted values from object
[`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md),
return by
[`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)

## Usage

``` r
# S3 method for class 'missSBM_fit'
fitted(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)

- ...:

  additional parameters for S3 compatibility.

## Value

A matrix of estimated probabilities of connection

## Examples

``` r
## a fit of the French political blogosphere, 30% of its dyads missing
adjMatrix <- missSBM::frenchblog2007 %>%
  igraph::as_adjacency_matrix(sparse = FALSE) %>%
  missSBM::observeNetwork(sampling = "dyad", parameters = 0.7)
collection <- estimateMissSBM(adjMatrix, 1:3, "dyad",
  control = missSBM_param(trace = FALSE, iterates = 0))
myModel <- collection$bestModel
fitted(myModel)[1:5, 1:5]
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,] 0.07674685 0.07674685 0.07674685 0.07674685 0.07674685
#> [2,] 0.07674685 0.07674685 0.07674685 0.07674685 0.07674685
#> [3,] 0.07674685 0.07674685 0.07674685 0.07674685 0.07674685
#> [4,] 0.07674685 0.07674685 0.07674685 0.07674685 0.07674685
#> [5,] 0.07674685 0.07674685 0.07674685 0.07674685 0.07674685
```
