# Changelog

## missSBM 1.2.0

This version changes the estimation of the missing-not-at-random designs
and the ICL: the results of `block-node`, `block-dyad`,
`double-standard` and `degree`, and the models selected for all the
designs, can differ from those of missSBM 1.1.0. To reproduce former
results, install missSBM 1.1.0 (e.g.,
`remotes::install_version("missSBM", "1.1.0")`).

### Changes in the estimation

- `block-node` and `block-dyad`: the sampling only depends on the
  blocks, so the missing dyads are integrated out given the blocks. The
  SBM part is fitted on the observed dyads, with the sampling term in
  the update of the blocks. This is a tighter bound of the same
  likelihood, and it removes a feedback of the imputations on the
  blocks, which could empty the blocks of unobserved nodes: the
  clustering is much more accurate when many dyads are missing.
- `double-standard`: the missing dyads are integrated out given the
  blocks too, and the connectivity and the sampling parameters (p0, p1)
  are estimated jointly. (p0, p1) are only identified through the
  differences of connectivity between the blocks: with a single block,
  the model is fitted with p0 = p1, i.e. as `dyad` sampling, and
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
  warns.
- `degree`: the variational EM maximizes the bound of the model again
  (Tabouy, Barbillon & Chiquet, 2020). It fits the SBM part on the
  observed dyads, then the sampling part, then all the parameters
  jointly, the blocks taking the imputed dyads into account. In the
  previous versions, the clustering ignored the sampling.
- The variational EM monitors and stops on the bound it maximizes
  (`loglik`), for all the designs.

### ICL

- The ICL is built on the observed data, the observed dyads and the
  sampling matrix: -2 times the ELBO at the MAP memberships plus the
  penalty. The missing dyads are integrated out exactly, except for
  `degree`, whose bound includes the imputed dyads and their entropy.
  This makes the ICL of different designs comparable, although some
  designs remain hard to tell apart (e.g. `block-node` and `degree` on
  an assortative network).
- The penalty counts the observed dyads for the connectivity instead of
  all the possible dyads.
- The ICL of the previous versions, computed on the network completed by
  the imputations, is no longer available: it rewarded the models that
  predict the missing dyads with certainty, and selected too many blocks
  when many dyads are missing.

### New fields

- `missSBM_fit$imputations`: the imputations of the missing dyads of
  `degree` (`NULL` for the other designs, which integrate them out;
  `imputedNetwork` gives the imputed network of any design).
- `fittedSBM$nbObservedDyads`; a column `stage` in the `monitoring` of a
  `degree` fit.

### Arguments and documentation

- A new vignette,
  [`vignette("case_study_frenchblog")`](https://grosssbm.github.io/missSBM/articles/case_study_frenchblog.md),
  follows a reduced form of the analysis of the French political
  blogosphere of the JSS article with the current version.
- [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
  and
  [`observeNetwork()`](https://grosssbm.github.io/missSBM/reference/observeNetwork.md)
  check their arguments with informative messages.
- [`?estimateMissSBM`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
  describes the parameters of every design and the model/design
  selection by the ICL.

### Memory and speed

- No n x n matrix is formed any more for the designs without covariates
  (network object, initial clustering, imputations, fit), and several
  computations run in C++: the memory used is several times smaller
  (building the network object takes about the size of its result,
  instead of five times the input), and the estimation is about twice as
  fast. Same results.

### Bug fixes

- `degree`: the imputations were not probabilities, the sampling
  parameters were estimated on wrong degrees (since 1.0.0), the
  variational parameter of the bound of Jaakkola and Jordan was
  truncated to \[0, 1\], which biased the sampling parameters (since
  2018), and the ICL used an approximation above the expected
  log-likelihood of the sampling instead of its bound.
- `block-dyad`: the expected log-likelihood of the sampling was the log
  of an expectation, which overestimated the objective and the ICL.
- `block-node` and `block-dyad`: the sampling parameters were not
  reordered with the blocks of the SBM, so that
  `fittedSampling$parameters` and
  [`coef()`](https://rdrr.io/r/stats/coef.html) did not match
  `fittedSBM`.
- `imputedNetwork` now gives the imputations at the final parameters
  (they came from the last iteration of the variational EM).
- Covariates in the SBM: the M-step could fail when two blocks had no
  edge between them, as the connectivity reached 0 exactly; it is now
  bounded away from 0 and 1, as without covariates.
- [`summary()`](https://rdrr.io/r/base/summary.html) of a
  `missSBM_collection` failed; it now prints the collection, as for a
  fit.
- The `war` data sets were stored by an old version of igraph, which
  printed a conversion message when used.
- Tests: the reference tests against `sbm` run sequentially
  (`nbCores = 1`), avoiding a segfault on macOS
  (`r-oldrel-macos-arm64`).

## missSBM 1.1.0 (2026-07-21)

CRAN release: 2026-07-21

### Major changes

- **Breaking**:
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)’s
  `control` must now be built with the new
  [`missSBM_param()`](https://grosssbm.github.io/missSBM/reference/missSBM_param.md)
  (one named, defaulted argument per option) instead of a raw
  `list(...)`, which now errors with a message pointing at the
  replacement.
- `missSBM_fit` now exposes
  [`split()`](https://rdrr.io/r/base/split.html),
  [`merge()`](https://rdrr.io/r/base/merge.html), `candidates_split()`
  and `candidates_merge()` as instance methods, previously inlined in
  `missSBM_collection`’s exploration logic; same search algorithm, now
  independently testable.
- new `polish(control)` (`missSBM_fit` and `missSBM_collection`):
  node-swap refinement after VEM convergence, fixing individually
  misclassified nodes that
  [`split()`](https://rdrr.io/r/base/split.html)/[`merge()`](https://rdrr.io/r/base/merge.html)
  cannot reach.
- requesting more blocks than a network actually supports can make VEM
  collapse one or more classes; this used to be silent, and split/merge
  exploration’s own repair of it could silently corrupt `vBlocks`’s
  bookkeeping (duplicated/missing entries, non-smooth ICL/ELBO in
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)). Now fixed
  and made visible: `missSBM_fit$repair(control)` recovers a collapsed
  fit ; new `occupiedBlocks`/ `degenerate` fields expose any remaining
  collapse; `bestModel` skips degenerate models when possible and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) marks them
  with a distinct point shape;
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)’s
  new `stopOnDegenerate`/`maxConsecutiveDegenerate` controls (default
  `TRUE`/2) stop forward exploration from growing further into a
  persistently unsupported range.
- new `estimate_chain()` (`missSBM_collection`) / `warmChain` control
  (default `FALSE`, opt-in): initializes each model by splitting the
  already-converged, smaller neighbor instead of an independent cold
  clustering, substantially reducing collapse in practice, at the cost
  of fitting `vBlocks` sequentially rather than in parallel.
- replace the NLopt/CCSAQ optimizer for the covariate connectivity
  parameters with a builtin Newton-Raphson solver (the M-step objective
  is concave, so it converges reliably in a handful of iterations);
  `nloptr` is no longer a dependency
- fix `partlyObservedNetwork$imputation()`: the fill value for missing
  dyads was biased low by including not-yet-imputed entries in its own
  computation. Known side effect: `doubleStandardSampling_fit`’s initial
  psi bootstrap can degenerate when the fill value exactly matches the
  empirical observed edge rate (refined away by the VEM loop that
  follows)

### Minor changes

- profiled
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
  and fixed several hot-path inefficiencies: avoidable `dense * sparse`
  promotion overhead, a redundant recomputation of the dense
  imputed-network matrix, and an
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) evaluating both of
  its branches unconditionally. Roughly 3x faster on our benchmark,
  bit-identical results
- [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
  now sorts/de-duplicates `vBlocks` if needed (with a warning), since
  exploration and chaining both assume it is strictly increasing.
- [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)’s
  `control` drops the `exploration` field: nobody used `"forward"`/
  `"backward"` alone at that level, only `"both"`/`"none"`, so
  `iterates` now doubles as the on/off switch (`0` disables exploration
  entirely). `missSBM_collection$explore()` no longer reads a stored
  `exploration` control field either; its `direction` argument (default
  `"both"`) is now a plain per-call parameter, still available for
  advanced/internal use.
- cap the number of merge candidates tried during backward exploration
  (`control$maxMergeCandidates`, default 30) instead of always trying
  every pair
- fix a crash and a closed-form algebra bug in the “degree” sampling
  design; parameter recovery for this design can still be biased under
  heavy missingness (known limitation)
- fix a consistency bug in `missSBM_fit$doVEM()`’s step-back: only the
  SBM was restored, not the sampling model or the current imputation
- speed up `getCovarArray()` and `kmeans_missSBM()`’s seeding
- remove unused `src/utils.h`
- `missSBM_fit`’s `polish()`, `repair()`, `candidates_split()` and
  `candidates_merge()` gain the same default `control` as `doVEM()`, so
  they can be called standalone without building a full control list
  first; `missSBM_collection`’s calls (which always pass an explicit
  control list) are unaffected ; `missSBM_collection`’s
  `estimate()`/`polish()`/`explore()` now share a control list stored at
  construction, no longer requiring `control` on every call.
- Factored `missSBM_collection`’s
  `explore_forward()`/`explore_backward()` into a shared private
  `explore_direction()`, and `missSBM_fit`’s
  `candidates_split()`/`candidates_merge()`’s trial-fit loop into a
  shared private `trial_fit_candidates()`; verified bit-identical output
  on the same fixed-seed scenario used throughout this changelog to
  diagnose VEM collapse

## missSBM 1.0.5 (2025-03-12)

CRAN release: 2025-03-13

- minor fix to comply with nlopt version 2.9.1 (NOCEDAL)

## missSBM 1.0.4 (2023-10-19)

CRAN release: 2023-10-24

- minor fix in package documentation due to evolution of roxygen2 7.0.0
  <https://github.com/r-lib/roxygen2/issues/1491>.

## missSBM 1.0.3

CRAN release: 2022-08-23

- Tiny adaptation due to new Matrix version 1.4-2
- fix for HTML5 in documentation, plus various typos

## missSBM 1.0.2

CRAN release: 2022-02-01

- Fix linking problem with new version of nloptR (2.0.0)
- Reference the JSS paper

## missSBM 1.0.1 \[2021-06-04\]

CRAN release: 2021-06-04

- less conservative tests to avoid random failure in CRAN checks
- tiny improvements in partlyObservedNetwork (less storage)

## missSBM 1.0.0 \[2021-05-25\]

CRAN release: 2021-05-25

- now rely on future_lapply for parallel computing (plan to be set by
  the user)
- faster model exploration (used to be called ‘smoothing’), now
  integrated by default in estimateMissSBM
- Use sparse Matrices to encode 0 and NAs
- complete rewriting of optimization routines (E and M steps) with C++
  armadillo routines
- Better initialization and embedded C++ kmeans implementation
- important bug fix in MAR case
- bug fix in inference on covariates
- bug fixed in blockDyad-sampling
- missSBM::SimpleSBM_fit_missSBM now inherits from from sbm::SimpleSBM
  rather than sbm::SimpleSBM_fit
- change field ‘$`netMatrix' to '`$networkData’ to comply with new
  interface in sbm
- defunct functions estimate, sample and simulate are no longer
  supported

## missSBM 0.3.0 \[2020-11-18\]

CRAN release: 2020-11-18

- changing interface after suggestion from JSS reviewers
- updated documentation
- interfacing with package sbm
  - change estimate to estimateMissSBM
  - change sample to observedNetwork
  - use sbm::sampleSimpleSBM instead of missSBM::simulate
  - export less R6 classes for simplification (internal use only)
- some bug fixes
- updated maintainer (<julien.chiquet@inra.fr> -\>
  <julien.chiquet@inrae.fr>)

## missSBM 0.2.2 \[2020-09-30\]

- unexporting sampledNetwork, only use internally
- merging prepare_data with estimate
- enhanced documentation
- moving ownership to großBM

## missSBM 0.2.1 \[2019-09-16\]

CRAN release: 2019-09-16

- added S3 methods for missSBM_fit, SBM_fit

## missSBM 0.2.0 \[2019-06-06\]

CRAN release: 2019-06-08

### significant changes:

- decent vignette
- faster tests
- many bugs corrected

## missSBM 0.1.0-9000 \[2019-02-26\]

- Added a `NEWS.md` file to track changes to the package.
