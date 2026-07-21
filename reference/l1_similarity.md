# L1-similarity

Compute l1-similarity between two vectors

## Usage

``` r
l1_similarity(x, y)
```

## Arguments

- x:

  a vector

- y:

  a vector

## Value

a vector equal to -abs(x-y)

## Examples

``` r
l1_similarity(1:5, 5:1)
#> [1] -4 -2  0 -2 -4
```
