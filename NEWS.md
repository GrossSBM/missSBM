# missSBM 1.1.0

## Major changes

- `missSBM_fit` now exposes `split()`, `merge()`, `candidates_split()` and `candidates_merge()`
  as instance methods, previously inlined in `missSBM_collection`'s exploration logic; same
  search algorithm, now independently testable.
- new `polish(control)` (`missSBM_fit` and `missSBM_collection`): node-swap refinement after VEM
  convergence, fixing individually misclassified nodes that `split()`/`merge()` cannot reach.
  Cheaper than split/merge exploration since it stays at a fixed number of blocks. Run
  automatically by `estimateMissSBM()` via its new `polish` control (default `TRUE`).
  `missSBM_collection`'s `estimate()`/`polish()`/`explore()` now share a control list stored at
  construction, no longer requiring `control` on every call.
- requesting more blocks than a network actually supports can make VEM collapse one or more
  classes; this used to be silent, and split/merge exploration's own repair of it could
  silently corrupt `vBlocks`'s bookkeeping (duplicated/missing entries, non-smooth ICL/ELBO in
  `plot()`). Now fixed and made visible: `missSBM_fit$repair(control)` recovers a collapsed fit
  (called automatically after every VEM fit and inside `polish()`); new `occupiedBlocks`/
  `degenerate` fields expose any remaining collapse; `bestModel` skips degenerate models when
  possible and `plot()` marks them with a distinct point shape; `estimateMissSBM()`'s new
  `stopOnDegenerate`/`maxConsecutiveDegenerate` controls (default `TRUE`/2) stop forward
  exploration from growing further into a persistently unsupported range.
- new `estimate_chain()` (`missSBM_collection`) / `warmChain` control (default `FALSE`, opt-in):
  initializes each model by splitting the already-converged, smaller neighbor instead of an
  independent cold clustering, substantially reducing collapse in practice, at the cost of
  fitting `vBlocks` sequentially rather than in parallel.
- `estimateMissSBM()` now sorts/de-duplicates `vBlocks` if needed (with a warning), since
  exploration and chaining both assume it is strictly increasing.
- replace the NLopt/CCSAQ optimizer for the covariate connectivity parameters with a builtin
  Newton-Raphson solver (the M-step objective is concave, so it converges reliably in a handful
  of iterations); `nloptr` is no longer a dependency
- profiled `estimateMissSBM()` and fixed several hot-path inefficiencies: avoidable
  `dense * sparse` promotion overhead, a redundant recomputation of the dense imputed-network
  matrix, and an `ifelse()` evaluating both of its branches unconditionally. Roughly 3x faster
  on our benchmark, bit-identical results
- fix `partlyObservedNetwork$imputation()`: the fill value for missing dyads was biased low by
  including not-yet-imputed entries in its own computation. Known side effect:
  `doubleStandardSampling_fit`'s initial psi bootstrap can degenerate when the fill value exactly
  matches the empirical observed edge rate (refined away by the VEM loop that follows)

## Minor changes

- cap the number of merge candidates tried during backward exploration
  (`control$maxMergeCandidates`, default 30) instead of always trying every pair
- fix a crash and a closed-form algebra bug in the "degree" sampling design; parameter recovery
  for this design can still be biased under heavy missingness (known limitation)
- fix a consistency bug in `missSBM_fit$doVEM()`'s step-back: only the SBM was restored, not the
  sampling model or the current imputation
- speed up `getCovarArray()` and `kmeans_missSBM()`'s seeding
- remove unused `src/utils.h`

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

