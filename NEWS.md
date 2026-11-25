# missSBM 1.2.0

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
  bound maximized by the algorithm, evaluated at the MAP memberships (each node in its most
  probable block, as the classification likelihood of Biernacki, Celeux & Govaert), plus the
  penalty; the same for the ICL of the fitted SBM alone (`fittedSBM$ICL`). For the designs
  whose sampling only depends on the blocks (`dyad`, `node`, `covar-*`, `block-node`,
  `block-dyad`), and for `double-standard`, the missing dyads are integrated out exactly. For
  `degree`, the bound includes the expected log-likelihood and the entropy of the imputations.
  (During the development of this version, the ICL was the bound minus the entropy of the soft
  memberships: under `degree`, whose imputations, computed from the memberships, make them
  near-hard, this rewarded `degree`, chosen in 20% to 75% of the cases under `node` and
  `block-node` sampling; at the MAP memberships, the designs are on an equal footing, #55.) The
  ICL of the previous
  versions added the expected log-likelihood of the imputed dyads, without their entropy, as if
  the missing dyads were classified: as the imputations are computed by the model, this
  rewards purer blocks with no data behind them, and favors more blocks, the more so as more
  dyads are missing. In simulations (`inst/simus/PR48-block-designs/`), the true number of
  blocks is selected in 95% to 100% of the cases in 11 settings out of 12 (it was 0% to 100%),
  and the ICL now picks the design that generated the data among `node` / `block-node` and
  `dyad` / `block-dyad` (it preferred `block-node` and `double-standard` in most cases). The
  selected models of all the designs can change.
- The penalty of the ICL built on the observed data counts the observed dyads for the
  connectivity of the SBM (and the effects of the covariates), estimated on them, instead of all
  the possible dyads (#51), as sbm does for a network with missing values. The sampling
  parameters still count all the dyads (or all the nodes): the sampling matrix is fully observed.
  The penalty of each connectivity parameter decreases by -log(rate of observed dyads), e.g. 1.9
  at 15% of observed dyads. In simulations (`inst/simus/PR51-penalty/`, 800 networks near the
  detectability threshold, Erdos-Renyi included), the selected number of blocks changes for 13
  networks, all toward the true one, and no structure is found in Erdos-Renyi networks with
  either penalty; the campaign of #48 is unchanged. New field `nbObservedDyads` of the fitted
  SBM.
- The ICL of the previous versions, built on the network completed by the imputations, is no
  longer available (it was kept as `missSBM_param(ICL = "completed")` during the development of
  this version, and is removed in #51): it has no justification, and it would not reproduce
  former results anyway, as the estimation of all the MNAR designs has changed. To reproduce
  former results, install missSBM 1.1.0 (e.g., `remotes::install_version("missSBM", "1.1.0")`).
  With a single block, `block-node` and `node`, or `block-dyad` and `dyad`, still give the same
  ICL.
- The VEM monitors (and stops on) the bound it maximizes, `loglik`, for all the designs. For
  `dyad`, `node` and `covar-*`, it monitored the expected log-likelihood of the network
  completed by the imputations of the SBM plus their entropy, which is not the bound: the fits
  of these designs can stop at a different iteration (in simulations, same selected number of
  blocks, ICL differences without direction).
- `degree` sampling (#55): the VEM maximizes the bound of the model again (factorized
  q(Z) q(Y_miss), bound of Jaakkola and Jordan), as in Tabouy, Barbillon & Chiquet (2020). Since
  missSBM 1.0.0, its SBM part was fitted on the observed dyads, ignoring the sampling (the
  clustering was that of `node` sampling), and its sampling part given the SBM part: no single
  objective was maximized. It is now fitted in three stages, each stopping on the objective it
  increases: the SBM part on the observed dyads, the sampling part given the SBM part, then the
  coordinate ascent of the bound in all the variables, the blocks and the connectivity being
  updated on the network completed by the imputations (new internal class
  `SimpleSBM_fit_degree`). The sampling part, which only involves the missing dyads, is iterated
  to convergence given the SBM part (in stage 2 and at each M-step of stage 3), as the EM of
  `double-standard`: with the default threshold, one step per iteration stopped it early. Started from the usual initialization, this coordinate ascent reaches
  worse local optima; from the first two stages, it increases the bound in every simulated
  network and keeps the partition. In simulations (`inst/simus/issue55-degree/`, n = 200), the
  ICL chooses `degree` when it generated the data in 100% of the cases for a medium or strong
  effect of the degree (0% to 20% before), 20% for a weak one, and in 5% of the cases without
  effect of the degree; in the campaign of #48, `degree` is chosen in 90% of the cases under
  `degree` sampling (0% before), the true design in 95% under `node` and 60% under `block-node`
  (the two designs are close for an assortative network). The `monitoring` of a `degree` fit
  has a column `stage`.
- `double-standard` sampling (#52): the missing dyads are also integrated out given the blocks,
  instead of being approximated by imputations independent of the blocks. A dyad between blocks
  q and l is missing with probability m_ql = pi_ql (1 - p1) + (1 - pi_ql) (1 - p0): this term
  enters the update of the blocks, and the connectivity and (p0, p1) are estimated jointly, by an
  EM on the missing dyads iterated to convergence at each M-step, in Q x Q. The previous version
  made one step of this EM per iteration of the VEM, which stopped long before convergence
  (after 2 iterations with the default threshold): p1 was underestimated (0.53 for 0.8, n = 200,
  Q = 3). In simulations (n = 200, Q = 3, 6 replicates per setting), the ARI goes from 0.73 to
  0.92, from 0.85 to 0.98 and from 0.96 to 0.99 for (p0, p1) = (0.6, 0.3), (0.5, 0.5),
  (0.2, 0.8), (p0, p1) are estimated without bias, and the true number of blocks is always
  selected. Among `dyad`, `block-dyad` and `double-standard`, the ICL now picks
  `double-standard` in 100% of the cases when it generated the data (60% before, see
  `inst/simus/PR48-block-designs/`). No imputation of each missing dyad is stored any more:
  `imputedNetwork` averages the posterior given the blocks over the memberships.
- `double-standard` with a single block: (p0, p1) are only identified through the differences of
  connectivity between the blocks. With a single block, the observed data only inform pi p1 and
  (1 - pi) p0: the model is fitted with p0 = p1, i.e. as dyad sampling (same ICL), with one
  sampling parameter in the penalty, and `estimateMissSBM()` warns when `vBlocks` contains 1.

## Memory and speed (#49)

- No n x n matrix is formed any more for the designs without covariates: the imputations are not
  computed along the VEM of the designs that integrate the missing dyads out (they are computed
  on demand by `imputedNetwork`, at the final parameters), and at the missing dyads only, in C++;
  the network object is built by reading the dense input in C++, without masks; the initial
  clustering uses no matrix of distances. At n = 4000 (`inst/simus/issue49-memory/`), building
  the network object takes +95 Mb instead of 5 times the input (+574 to 715 Mb), the initial
  clustering +91 to 128 Mb instead of +526 to 672 Mb, and the fit of `node`, `block-dyad` and
  `double-standard` +4 to 110 Mb instead of +301 to 637 Mb; the whole estimation is about twice
  as fast. Same results (up to rounding for the imputations).
- `degree` (#57): the update of its imputations runs in C++, in one pass over the missing dyads
  (it allocated about ten vectors of their size at each of the many iterations of its sampling
  part), and the E-step no longer transposes the sparse matrices. At n = 4000, its fit goes from
  6.7 s to 4.2 s (installed package). Same results, up to rounding.

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
- `degree` sampling, bound of Jaakkola and Jordan: the variational parameter zeta_i, at its
  optimum sqrt(E[(a + b D_i)^2]), often above 1, was truncated to [eps, 1 - eps] like a
  probability (since 2018). The estimates of (a, b) and the imputations were biased, the more so
  as the effect of the degree is strong (b underestimated by 15% for a strong effect, n = 200; now
  0.6%; `inst/simus/issue55-degree/`). Also, the expected log-likelihood of the sampling, used in
  the ICL, plugged the expected degrees into log logistic(a + b D_i), which is above its
  expectation (Jensen), instead of the bound of Jaakkola and Jordan that the updates maximize: the
  ICL of `degree` now uses this bound.
- `block-dyad` sampling: the expectation of the log-likelihood of the sampling was computed as
  the log of the expected sampling probabilities, log(Z psi Z'), instead of the expectation of
  the log (since 2018). The updates of the parameters were right, but the objective monitored by
  the VEM and the ICL were overestimated (by 5 to 20 on the ICL of networks with 200 nodes).
  The estimates are unchanged in our simulations, but the selected number of blocks can change.
- `block-node` and `block-dyad` sampling: the sampling parameters, indexed by the blocks, were
  not permuted with the blocks of the SBM when these are reordered at the end of the VEM, so that
  `fittedSampling$parameters` (and `coef()`) did not match the blocks of `fittedSBM`, and
  `polish()` added a permuted sampling term to its gains. The ICL was not affected. A side effect:
  this permuted term acted as a random perturbation that could get `polish()` out of a local
  optimum, so that a fit at a fixed number of blocks without exploration can now end elsewhere.

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
