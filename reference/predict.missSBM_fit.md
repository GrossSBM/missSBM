# Prediction of a [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md) (i.e. network with imputed missing dyads)

Prediction of a
[`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)
(i.e. network with imputed missing dyads)

## Usage

``` r
# S3 method for class 'missSBM_fit'
predict(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)

- ...:

  additional parameters for S3 compatibility.

## Value

an adjacency matrix between pairs of nodes. Missing dyads are imputed
with their expected values, i.e. by there estimated probabilities of
connection under the missing SBM.

## Examples

``` r
## a fit of the French political blogosphere, 30% of its dyads missing
adjMatrix <- missSBM::frenchblog2007 %>%
  igraph::as_adjacency_matrix(sparse = FALSE) %>%
  missSBM::observeNetwork(sampling = "dyad", parameters = 0.7)
collection <- estimateMissSBM(adjMatrix, 1:3, "dyad",
  control = missSBM_param(trace = FALSE, iterates = 0))
myModel <- collection$bestModel
## the observed dyads, and the imputations of the missing ones
predict(myModel)[1:5, 1:5]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                                                            
#> [1,] .          1.00000000 0.07618759 0.07618759 0.07618759
#> [2,] 1.00000000 .          1.00000000 0.07618759 1.00000000
#> [3,] 0.07618759 1.00000000 .          1.00000000 1.00000000
#> [4,] 0.07618759 0.07618759 1.00000000 .          1.00000000
#> [5,] 0.07618759 1.00000000 1.00000000 1.00000000 .         
```
