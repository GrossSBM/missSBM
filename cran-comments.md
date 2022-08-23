
This submission is in response to a request from a CRAN maintainer.
The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN.

## missSBM 1.0.3

- Adaptation due to new Matrix version 1.4-2
- fix for HTML5 in documentation

## Tested environments

Success

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

- tested remotely with win-builder (R version 4.0.5, R version 4.1.0,  R unstable)
- tested remotely with R-hub 
     * R-release: Ubuntu, Solaris, Windows Server 2008 R2 SP1, macOS HighSierra

Failure: Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results

On Ubuntu 20.04, Debian Gcc and Fedora Linux 1 NOTE about installation size

* checking installed package size ... NOTE
  installed size is  7.2Mb
  sub-directories of 1Mb or more:
    libs   5.1Mb

