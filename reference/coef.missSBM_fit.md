# Extract model coefficients

Extracts model coefficients from objects
[`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)
returned by
[`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)

## Usage

``` r
# S3 method for class 'missSBM_fit'
coef(
  object,
  type = c("mixture", "connectivity", "covariates", "sampling"),
  ...
)
```

## Arguments

- object:

  an R6 object with class
  [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)

- type:

  type of parameter that should be extracted. Either "mixture"
  (default), "connectivity", "covariates" or "sampling"

- ...:

  additional parameters for S3 compatibility. Not used

## Value

A vector or matrix of coefficients extracted from the missSBM_fit model.

## Examples

``` r
## a fit of the French political blogosphere, 30% of its dyads missing
adjMatrix <- missSBM::frenchblog2007 %>%
  igraph::as_adjacency_matrix(sparse = FALSE) %>%
  missSBM::observeNetwork(sampling = "dyad", parameters = 0.7)
collection <- estimateMissSBM(adjMatrix, 1:3, "dyad",
  control = missSBM_param(trace = FALSE, iterates = 0))
myModel <- collection$bestModel
coef(myModel, "mixture")
#> [1] 0.64676714 0.01528746 0.33794539
coef(myModel, "connectivity")
#> $mean
#>              [,1]         [,2]         [,3]
#> [1,] 1.218759e-01 2.220446e-16 1.293078e-02
#> [2,] 2.220446e-16 1.000000e+00 2.220446e-16
#> [3,] 1.293078e-02 2.220446e-16 1.600133e-01
#> 
coef(myModel, "sampling")
#> [1] 0.7009419
```
