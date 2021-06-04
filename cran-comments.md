
## missSBM 1.0.1

Minor version with less conservative unitary tests to avoid random error during (CRAN) checks

## Tested environments

Rhub
Windows R Under development (unstable)
Windows R version 4.0.5
Ubuntu Linux 20.04.1 LTS, R-release, GCC
Fedora Linux, R-devel, clang, gfortran
	
- tested remotely with github action (macOS, Ubuntu, Windows release)
- tested remotely with win-builder (OK on release, old, failure on R-devel)
- tested remotely with R-hub (release Ubuntu, devel Fedora, release Solaris, release macOS High Sierra)
- tested locally on Ubuntu 20.04, R 4.1.0

## R CMD check results

On Ubuntu 20.04, Fedora Linux 1 NOTE about installed size

── R CMD check results ────────────────────────────────────── missSBM 1.0.0 ────

* checking installed package size ... NOTE
  installed size is 15.7Mb
  sub-directories of 1Mb or more:
    libs    8.7Mb

0 errors ✓ | 0 warnings ✓ | 1 note 

R CMD check succeeded
