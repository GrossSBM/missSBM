
# missSBM 1.0.4 (2023-10-24)

minor fix in package documentation due to evolution of roxygen2 7.0.0 <https://github.com/r-lib/roxygen2/issues/1491>, as asked by CRAN maintainer.

## Tested environments

* tested locally on Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with R-hub 
  - Windows Server 2022, R-devel, 64 bit
  - Fedora Linux, R-devel, clang, gfortran
	- Ubuntu Linux 20.04.1 LTS, R-release, GCC

* tested remotely with github-action
  - Linux ubuntu 22.04, R-release (github-action)
  - Linux ubuntu 22.04, R-oldrel (github-action)
  - Linux ubuntu 22.04, R-devel (github-action)
  - Windows Server 2022, R-release, 64 bit
  - macOS 12, R-release (github action)

- tested remotely with win-builder (R version 4.3.1, R version 4.2.3, R unstable)

## R CMD check results

1 NOTE about installation size

* checking installed package size ... NOTE
  installed size is  7.2Mb
  sub-directories of 1Mb or more:
    libs   5.1Mb
