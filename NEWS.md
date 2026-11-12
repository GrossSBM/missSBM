# missSBM (development version)

## Changes in the estimation

- `block-node` and `block-dyad` sampling: the SBM part is now fitted on the observed dyads,
  with the sampling term in the update of the blocks, as was already the case for `node` and
  `dyad`. Under these designs the sampling only depends on the blocks, so that the missing
  dyads are ignorable given the blocks: they are now integrated out exactly, instead of being
  approximated by a variational distribution independent of the blocks, q(Z) q(Y_miss). This is
  a tighter lower bound of the same likelihood, not a new criterion. It also removes a feedback
  loop: the imputations of the missing dyads, computed from the blocks, entered the update of
  the blocks, and the dyads between two unobserved nodes, never observed, fed back on themselves
  (the connectivity of a block of unobserved nodes could drift to 0, and classes collapsed when
  exploring larger numbers of blocks, e.g., the 196 blogs of `frenchblog2007` up to 18 blocks).
  The results of these two designs change. In simulations (n = 200, Q = 3), the mean ARI goes
  from 0.65 to 0.94 for `block-node` (42% of observed nodes) and from 0.32 to 0.88 for
  `block-dyad` (30% of observed dyads).
- The ICL is now built on the observed data, for all the designs: the observed dyads and the
  sampling matrix, the missing dyads being integrated out. It is -2 times the variational lower
  bound maximized by the algorithm, without the entropy of the blocks, plus the penalty: the
  same relation between the ICL and the bound as for a fully observed SBM. For the designs
  whose sampling only depends on the blocks (`dyad`, `node`, `covar-*`, `block-node`,
  `block-dyad`), the missing dyads are integrated out exactly. For `degree` and
  `double-standard`, the bound includes the entropy of the imputations. The ICL of the previous
  versions added the expected log-likelihood of the imputed dyads, without their entropy, as if
  the missing dyads were classified: as the imputations are computed by the model, this
  rewards purer blocks with no data behind them, and favors more blocks, the more so as more
  dyads are missing. In simulations (`inst/simus/PR48-block-designs/`), the true number of
  blocks is selected in 95% to 100% of the cases in 11 settings out of 12 (it was 0% to 100%),
  and the ICL now picks the design that generated the data among `node` / `block-node` and
  `dyad` / `block-dyad` (it preferred `block-node` and `double-standard` in most cases). The
  selected models of all the designs can change.
- `missSBM_param(ICL = "completed")` restores the ICL of the previous versions, to reproduce
  former results. Both are available on each model (fields `ICLobserved` and `ICLcompleted` of
  `missSBM_fit`). With a single block, `block-node` and `node`, or `block-dyad` and `dyad`,
  still give the same ICL, with both.
- Only `double-standard`, whose sampling depends on the value of each dyad, still fits the SBM
  part on the imputations of the missing dyads (see #52).

## Bug fixes

- `degree` sampling: the imputations of the missing dyads were not probabilities (the
  logistic function had been dropped from their update in 2018, so that they were the
  linear predictor truncated to [eps, 1 - eps]), and the sampling parameters were estimated
  on wrong degrees (the observed edges were left out of the imputed network since 1.0.0, and
  only one triangle of an undirected network was read). The estimates of the sampling
  parameters and the imputations of this design are now correct; the SBM part is still fitted
  on the observed dyads, which is more robust than plugging the imputations of the dyads
  between unobserved nodes, as these are never observed.
- `imputedNetwork` no longer doubles the imputations with the `degree` sampling.
- `block-dyad` sampling: the expectation of the log-likelihood of the sampling was computed as
  the log of the expected sampling probabilities, log(Z psi Z'), instead of the expectation of
  the log (since 2018). The updates of the parameters were right, but the objective monitored by
  the VEM and the ICL were overestimated (by 5 to 20 on the ICL of networks with 200 nodes).
  The estimates are unchanged in our simulations, but the selected number of blocks can change.

## Minor changes

- `degree` sampling: the sampling part is computed on the sparse pattern of the missing dyads,
  in O(n + number of missing dyads) memory instead of several dense n x n matrices. Same results
  up to rounding errors (checked against the dense formulas in the tests).
- The observed nodes of an undirected network are found without building `R | t(R)`, a nearly
  full sparse matrix when most dyads are observed. This was done for every fit with a
  node-centered sampling, including each candidate of the exploration: about 18% faster overall
  on a network with 2000 nodes, same results.

# missSBM 1.1.1

## Bug fixes

- Force sequential execution (`nbCores = 1`) in the `blockmodels`-based
  reference tests, avoiding a segfault caused by `blockmodels`'s internal
  `parallel::mclapply` use, which is unsafe when combined with macOS's
  Accelerate/vecLib BLAS-LAPACK framework (observed on CRAN's
  r-oldrel-macos-arm64 check machine).

# missSBM 1.1.0 (2026-07-21)

## Major changes

- **Breaking**: `estimateMissSBM()`'s `control` must now be built with the new `missSBM_param()`
  (one named, defaulted argument per option) instead of a raw `list(...)`, which
  now errors with a message pointing at the replacement.
- `missSBM_fit` now exposes `split()`, `merge()`, `candidates_split()` and `candidates_merge()`
  as instance methods, previously inlined in `missSBM_collection`'s exploration logic; same
  search algorithm, now independently testable.
- new `polish(control)` (`missSBM_fit` and `missSBM_collection`): node-swap refinement after VEM
  convergence, fixing individually misclassified nodes that `split()`/`merge()` cannot reach.
- requesting more blocks than a network actually supports can make VEM collapse one or more
  classes; this used to be silent, and split/merge exploration's own repair of it could
  silently corrupt `vBlocks`'s bookkeeping (duplicated/missing entries, non-smooth ICL/ELBO in
  `plot()`). Now fixed and made visible: `missSBM_fit$repair(control)` recovers a collapsed fit ;
  new `occupiedBlocks`/ `degenerate` fields expose any remaining collapse;
  `bestModel` skips degenerate models when possible and `plot()` marks them with a distinct point
  shape; `estimateMissSBM()`'s new `stopOnDegenerate`/`maxConsecutiveDegenerate` controls
  (default `TRUE`/2) stop forward exploration from growing further into a persistently unsupported range.
- new `estimate_chain()` (`missSBM_collection`) / `warmChain` control (default `FALSE`, opt-in):
  initializes each model by splitting the already-converged, smaller neighbor instead of an
  independent cold clustering, substantially reducing collapse in practice, at the cost of
  fitting `vBlocks` sequentially rather than in parallel.
- replace the NLopt/CCSAQ optimizer for the covariate connectivity parameters with a builtin
  Newton-Raphson solver (the M-step objective is concave, so it converges reliably in a handful
  of iterations); `nloptr` is no longer a dependency
- fix `partlyObservedNetwork$imputation()`: the fill value for missing dyads was biased low by
  including not-yet-imputed entries in its own computation. Known side effect:
  `doubleStandardSampling_fit`'s initial psi bootstrap can degenerate when the fill value exactly
  matches the empirical observed edge rate (refined away by the VEM loop that follows)

## Minor changes

- profiled `estimateMissSBM()` and fixed several hot-path inefficiencies: avoidable
  `dense * sparse` promotion overhead, a redundant recomputation of the dense imputed-network
  matrix, and an `ifelse()` evaluating both of its branches unconditionally. Roughly 3x faster
  on our benchmark, bit-identical results
- `estimateMissSBM()` now sorts/de-duplicates `vBlocks` if needed (with a warning), since
  exploration and chaining both assume it is strictly increasing.
- `estimateMissSBM()`'s `control` drops the `exploration` field: nobody used `"forward"`/
  `"backward"` alone at that level, only `"both"`/`"none"`, so `iterates` now doubles as the
  on/off switch (`0` disables exploration entirely). `missSBM_collection$explore()` no longer
  reads a stored `exploration` control field either; its `direction` argument (default `"both"`)
  is now a plain per-call parameter, still available for advanced/internal use.
- cap the number of merge candidates tried during backward exploration
  (`control$maxMergeCandidates`, default 30) instead of always trying every pair
- fix a crash and a closed-form algebra bug in the "degree" sampling design; parameter recovery
  for this design can still be biased under heavy missingness (known limitation)
- fix a consistency bug in `missSBM_fit$doVEM()`'s step-back: only the SBM was restored, not the
  sampling model or the current imputation
- speed up `getCovarArray()` and `kmeans_missSBM()`'s seeding
- remove unused `src/utils.h`
- `missSBM_fit`'s `polish()`, `repair()`, `candidates_split()` and `candidates_merge()` gain the
  same default `control` as `doVEM()`, so they can be called standalone without building a full
  control list first; `missSBM_collection`'s calls (which always pass an explicit control list)
  are unaffected ; `missSBM_collection`'s `estimate()`/`polish()`/`explore()` now share a control
  list stored at construction, no longer requiring `control` on every call.
- Factored `missSBM_collection`'s `explore_forward()`/`explore_backward()` into a shared private
  `explore_direction()`, and `missSBM_fit`'s `candidates_split()`/`candidates_merge()`'s trial-fit
  loop into a shared private `trial_fit_candidates()`; verified bit-identical output on the same
  fixed-seed scenario used throughout this changelog to diagnose VEM collapse

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
