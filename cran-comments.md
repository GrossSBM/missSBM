
This submission is in response to a request from a CRAN maintainer.

## missSBM 1.0.2

Fix linking problem with new version of nloptR

## Tested environments

Success

- tested remotely with github action (R release on macOS 10.15, Ubuntu 20.04, Debian gcc, Windows Server 2019)
- tested remotely with win-builder (R version 4.0.5, R version 4.1.0,  R unstable)
- tested remotely with R-hub 
     * R-release: Ubuntu, Solaris, Windows Server 2008 R2 SP1, macOS HighSierra
     * R devel: Fedora,
     * Windows Server 2008 R2 SP1 R old-release
- tested locally on Ubuntu 20.04, R 4.1.0

Failure: Windows Server 2008 R2 SP1, R-devel, 32/64 bit

## R CMD check results

On Ubuntu 20.04, Debian Gcc and Fedora Linux 1 NOTE about installation size

* checking installed package size ... NOTE
  installed size is 10.8Mb
  sub-directories of 1Mb or more:
    libs   8.7Mb

