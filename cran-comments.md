

# missSBM 1.0.4 (2023-10-19)

minor fix in package documentation due to evolution of roxygen2 7.0.0 <https://github.com/r-lib/roxygen2/issues/1491>, as asked by CRAN maintainer

## Tested environments

* tested locally on Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with R-hub 
  - Windows Server 2022, R-devel, 64 bit
  - Fedora Linux, R-devel, clang, gfortran
	- Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with github-action
  - Linux ubuntu 20.04, R-release (github-action)
  - Linux ubuntu 20.04, R-oldrel (github-action)
  - Linux ubuntu 20.04, R-devel (github-action)
  - Windows Server 2022, R-release, 64 bit
  - macOS Big Sur 11, R-release (github action)

- tested remotely with win-builder (R version 4.1.2, R version 4.1.3)

## R CMD check results

On Ubuntu and Fedora Linux 1 NOTE about installation size

* checking installed package size ... NOTE
  installed size is  7.2Mb
  sub-directories of 1Mb or more:
    libs   5.1Mb
