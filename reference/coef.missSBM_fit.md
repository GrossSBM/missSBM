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
