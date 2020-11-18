
## missSBM 0.3.0

  - changing interface after suggestion from JSS reviewers
  - updated documentation
  - interfacing with package sbm
    - change estimate to estimateMissSBM
    - change sample to observedNetwork
    - use sbm::sampleSimpleSBM instead of missSBM::simulate
    - export less R6 classes for simplification (internal use only)
  - some bug fixes
  - updated maintainer (julien.chiquet@inra.fr -> julien.chiquet@inrae.fr)

## Tested environments

- tested remotely with github action (macOS, ubuntu, Windows release, unbuntu devel)
- tested remotely with win-builder (devel)
- tested remotely with R-hub (release Ubuntu, Fedora)
- tested locally on Ubuntu 20.04, R 4.0.2

## R CMD check results

── R CMD check results ────────────────────────────────────── missSBM 0.3.1 ────
Duration: 1m 15.5s

> checking installed package size ... NOTE
    installed size is  6.2Mb
    sub-directories of 1Mb or more:
      libs   4.1Mb

0 errors ✓ | 0 warnings ✓ | 1 note x

R CMD check succeeded
