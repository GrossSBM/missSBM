
# missSBM 1.0.5 (2025-03-12)

minor fix to comply with nlopt version 2.9.1 (NOCEDAL)

## Tested environments

* tested locally on Ubuntu Linux 24.04.2 LTS, R-release, GCC

* tested remotely with win-builder 
  - Windows Server 2022, R 4.3.3, 64 bit
  - Windows Server 2022, R 4.4.3, 64 bit
  - Windows Server 2022, R-devel, 64 bit

* tested remotely with github-action
  - Linux Ubuntu 24.04, R-release
  - Linux Ubuntu 24.04, R-devel 
  - Windows Server 2022, R-release, 64 bit
  - Windows Server 2022, R-oldrel , 64 bit
  - macOS Sonoma 14, R-release 
  - macOS Sonoma 14, R-oldrel 

## R CMD check results

1 NOTE about installation size

N  checking installed package size ...
     installed size is  8.7Mb
     sub-directories of 1Mb or more:
       libs   6.6Mb
