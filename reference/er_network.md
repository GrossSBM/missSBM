# ER ego centered network

A dataset containing the weighted PPI network centered around the ESR1
(ER) protein

## Usage

``` r
er_network
```

## Format

A sparse symmetric matrix with 741 rows and 741 columns `ESR1`

## Source

<https://string-db.org/>

## Examples

``` r
data("er_network")
class(er_network)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"
```
