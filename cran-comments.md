
## missSBM 1.0.0

  - compatibility with the new CRAN version of packge 'sbm' 
  - now rely on future_lapply  for parallel computing (plan to be set by the user)
  - faster model exploration (used to be called 'smoothing'), now integrated by default in estimateMissSBM
  - Use sparse Matrices to encode 0 and NAs
  - complete rewriting of optimization routines (E and M steps) with C++ armadillo routines
  - Better initialization and embedded C++ kmeans implementation
  - important bug fix in MAR case
  - bug fix in inference on covariates
  - bug fixed in blockDyad-sampling
  - missSBM::SimpleSBM_fit_missSBM now inherits from from sbm::SimpleSBM rather than sbm::SimpleSBM_fit
  - change field '$netMatrix' to '$networkData' to comply with new interface in sbm
  - defunct functions estimate, sample and simulate are no longer supported

## Tested environments

- tested remotely with github action (macOS, Ubuntu, Windows release)
- tested remotely with win-builder (release, devel, old)
- tested remotely with R-hub (release Ubuntu, Fedora, Solaris, mac/CRAN setup)

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
