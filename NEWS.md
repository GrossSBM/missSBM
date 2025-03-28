# missSBM 1.0.5 (2025-03-12)

- minor fix to comply with nlopt version 2.9.1 (NOCEDAL)

# missSBM 1.0.4 (2023-10-19)

- minor fix in package documentation due to evolution of roxygen2 7.0.0 <https://github.com/r-lib/roxygen2/issues/1491>.

# missSBM 1.0.3

- Tiny adaptation due to new Matrix version 1.4-2
- fix for HTML5 in documentation, plus various typos

# missSBM 1.0.2

  - Fix linking problem with new version of nloptR (2.0.0)
  - Reference the JSS paper
  
# missSBM 1.0.1 [2021-06-04]

  - less conservative tests to avoid random failure in CRAN checks
  - tiny improvements in partlyObservedNetwork (less storage)

# missSBM 1.0.0 [2021-05-25]

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

# missSBM 0.3.0 [2020-11-18]

  - changing interface after suggestion from JSS reviewers
  - updated documentation
  - interfacing with package sbm
    - change estimate to estimateMissSBM
    - change sample to observedNetwork
    - use sbm::sampleSimpleSBM instead of missSBM::simulate
    - export less R6 classes for simplification (internal use only)
  - some bug fixes
  - updated maintainer (julien.chiquet@inra.fr -> julien.chiquet@inrae.fr)

# missSBM 0.2.2 [2020-09-30]

  - unexporting sampledNetwork, only use internally
  - merging prepare_data with estimate
  - enhanced documentation
  - moving ownership to großBM

# missSBM 0.2.1 [2019-09-16]
 
  - added S3 methods for missSBM_fit, SBM_fit

# missSBM 0.2.0 [2019-06-06]

## significant changes:
  - decent vignette
  - faster tests
  - many bugs corrected

# missSBM 0.1.0-9000 [2019-02-26]

* Added a `NEWS.md` file to track changes to the package.

