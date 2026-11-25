# Summary method for a [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)

Summary method for a
[`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)

## Usage

``` r
# S3 method for class 'missSBM_fit'
summary(object, ...)
```

## Arguments

- object:

  an R6 object with class
  [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)

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
myModel <- collection$bestModel
summary(myModel)
#> missSBM-fit
#> ==================================================================
#> Structure for storing an SBM fitted under missing data condition  
#> ==================================================================
#> * Useful fields (first 2 special objects themselves with methods) 
#>   $fittedSBM (the adjusted stochastic block model)                
#>   $fittedSampling (the estimated sampling process)                
#>   $imputedNetwork (the adjacency matrix with imputed values)      
#>   $monitoring, $ICL, $loglik, $vExpec, $penalty                  
#> * S3 methods                                                      
#>   plot, coef, fitted, predict, print                              
```
