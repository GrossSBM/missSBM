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
