
This submission is in response to a request from a CRAN maintainer.
The DOI in the CITATION is for a new JSS publication that will be registered after publication on CRAN.

## missSBM 1.0.3

- Adaptation due to new Matrix version 1.4-2
- fix for HTML5 in documentation

## missSBM 1.0.2

- Fix linking problem with new version of nloptR (2.0.0)
- Reference to JSS paper

## Tested environments

Success

- Windows Server 2022, R-devel, 64 bit
- Windows Server 2019, R-release (github action)
- macOS 10.15, R-release (github action)
- Linux ubuntu 20.04, R-release (github-action)
- Linux ubuntu 20.04, R-devel (github-action)
- Linux ubuntu 20.04, R-oldrel (github-action)

- tested remotely with win-builder (R version 4.0.5, R version 4.1.0,  R unstable)
- tested remotely with R-hub 
     * R-release: Ubuntu, Solaris, Windows Server 2008 R2 SP1, macOS HighSierra
     * R devel: Fedora,
     * Windows Server 2008 R2 SP1 R old-release

Failure: Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results

On Ubuntu 20.04, Debian Gcc and Fedora Linux 1 NOTE about installation size

* checking installed package size ... NOTE
  installed size is 10.8Mb
  sub-directories of 1Mb or more:
    libs   8.7Mb

