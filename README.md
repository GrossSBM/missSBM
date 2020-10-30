
# missSBM: Handling missing data in Stochastic Block Models

<!-- badges: start -->

[![website](https://github.com/GrossSBM/missSBM/workflows/pkgdown/badge.svg)](https://grosssbm.github.io/missSBM/)
[![R-CMD-check](https://github.com/grosssbm/missSBM/workflows/R-CMD-check/badge.svg)](https://github.com/grosssbm/missSBM/actions)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/missSBM)](https://cran.r-project.org/package=missSBM)
[![](https://img.shields.io/github/last-commit/grossSBM/missSBM.svg)](https://github.com/GrossSBM/missSBM/commits/master)
[![Codecov test
coverage](https://codecov.io/gh/GrossSBM/missSBM/branch/master/graph/badge.svg)](https://codecov.io/gh/GrossSBM/missSBM?branch=master)

<!-- badges: end -->

> When a network is partially observed (here, NAs in the adjacency
> matrix rather than 1 or 0 due to missing information between node
> pairs), it is possible to account for the underlying process that
> generates those NAs. ‘missSBM’ adjusts the popular stochastic block
> model from network data observed under various missing data
> conditions, as described in Tabouy, Barbillon and Chiquet (2019)
> [10.1080/01621459.2018.1562934](https://doi.org/10.1080/01621459.2018.1562934).

## Installation

The Last CRAN version is available via

``` r
install.packages("missSBM")
```

The development version is available via

``` r
devtools::install_github("grossSBM/missSBM")
```

## Reference

Please cite our work using the following reference:

Timothée Tabouy, Pierre Barbillon & Julien Chiquet (2019) “Variational
Inference for Stochastic Block Models from Sampled Data”, *Journal of
the American Statistical Association*, DOI:
[10.1080/01621459.2018.1562934](https://doi.org/10.1080/01621459.2018.1562934)
