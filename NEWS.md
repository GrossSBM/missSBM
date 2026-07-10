# missSBM 1.0.6

- `missSBM_fit` now exposes `split()`, `merge()`, `candidates_split()` and
  `candidates_merge()` as instance methods, mirroring the architecture of a
  sibling project (normalblockr's `NormalBlockBase`): the split/merge search
  algorithm itself is unchanged (a spectral bipartition of a cluster's induced
  sub-network for splitting, direct relabeling for merging), only its code
  organization moved from being inlined in `missSBM_collection`'s
  `explore_forward()`/`explore_backward()` to living on the model itself,
  which is now independently testable and readable. Because the exact
  boundaries of parallel (`future_lapply`) evaluation shifted slightly, the
  precise RNG draws consumed during exploration differ from before (verified:
  same ICL at every q on a fixed-seed benchmark except where a split/merge
  candidate was actually applied, where it differs by a negligible amount,
  never worse)
- add SQUAREM acceleration (Varadhan & Roland, 2008) to the VEM of
  `SimpleSBM_fit_noCov` (the MAR, no-covariate case): every two plain VEM
  steps, attempts an extrapolated step in an unconstrained reparametrization
  of `(theta, pi)` (`logit`/`log`, so any extrapolated point maps back to a
  feasible probability -- no separate feasibility guard needed), stabilized
  by one more E-step/M-step, accepted only if it does at least as well as
  plain VEM. Measured 1.6x-2.7x fewer iterations on harder (larger Q, noisy
  initial clustering) synthetic fits, matching the same fixed point in most
  cases; being a genuinely different optimization trajectory on a non-convex
  objective, it can occasionally converge to a different (and occasionally
  slightly worse) local optimum than plain VEM would have -- an inherent,
  known characteristic of extrapolation-based EM acceleration on multimodal
  likelihoods, observed in about 1 run out of 8 in that same experiment.
  Ported from and validated the same way as a similar acceleration in a
  sibling project (normalblockr); other model variants (with covariates,
  MNAR sampling designs, `missSBM_fit`'s composite SBM+sampling model) are
  deliberately not accelerated yet, pending the same per-class validation
- cap the number of merge candidates tried during backward exploration
  (`control$maxMergeCandidates`, default 30) instead of always trying all
  `choose(q, 2)` pairs: beyond the cap, only the pairs with the most similar
  fitted connectivity profiles are tried (merging two very different blocks
  is rarely competitive anyway). Ported from a similar idea validated in a
  sibling project (normalblockr's `candidates_merge()`)
- fix `partlyObservedNetwork$imputation()`: the "average"/"median" fill value for missing
  dyads was computed over the whole adjacency matrix, including the not-yet-imputed
  (still-zero) missing entries themselves, biasing it low; it is now computed from the
  observed entries only, and filled in via a sparse addition instead of a large-vector
  index assignment (`Matrix::dgCMatrix` `[<-` on many indices is a known-slow pattern)
- as a side effect of the fix above, `doubleStandardSampling_fit`'s one-shot initial
  bootstrap of psi (built on `imputation("average")`) is now more exposed to a structural
  degeneracy: whenever the fill value used equals the exact empirical observed edge rate,
  psi[1] and psi[2] are mathematically forced to coincide. This only affects the initial,
  non-iterated estimate (refined by the VEM loop afterwards, unaffected in practice) --
  known limitation, not addressed in this release
- fix a crash in the "degree" sampling design (`degreeSampling_fit` referenced a
  non-existent field on `partlyObservedNetwork`)
- fix the closed-form update of the "degree" sampling parameters (wrong coefficients
  in the Jaakkola-Jordan variational M-step); note that parameter recovery for this
  design can still be biased under heavy missingness, this is a known limitation
- fix a consistency bug in `missSBM_fit$doVEM()`: when the variational EM stepped back
  after a decrease of the objective, only the fitted SBM was restored, not the fitted
  sampling model nor the current imputation, leaving them out of sync; both `doVEM()`
  implementations (`SimpleSBM_fit`, `missSBM_fit`) now share a common driver and use a
  lightweight state snapshot instead of a full clone of the (possibly large) SBM object
- speed up `getCovarArray()` (used to build the covariate similarity array) with a fully
  vectorized fast path for the default `l1_similarity`, replacing an O(N^2) R-level loop
- speed up `kmeans_missSBM()`'s farthest-point seeding: it now maintains a running
  per-point distance to the nearest already-chosen centroid instead of recomputing it
  from scratch at every iteration (O(k N) instead of O(k^2 N))
- replace the NLopt/CCSAQ optimizer used to fit the covariate connectivity parameters
  (`M_step_sparse_bernoulli_covariates`) with a builtin Newton-Raphson solver: the M-step
  objective is a weighted logistic regression, globally concave in `(Gamma, beta)`, so
  Newton converges reliably in a handful of iterations (typically 4-6, against 20-60 for
  CCSAQ at comparable precision) instead of relying on a generic black-box optimizer.
  The symmetric (undirected) case now optimizes the actual `Q(Q+1)/2` free parameters
  directly (instead of a redundant `Q x Q` matrix with a post-hoc symmetrized gradient),
  guaranteeing exact symmetry at every iteration. A Levenberg-Marquardt-style damping with
  backtracking line search keeps the solver robust when the Hessian is near-singular,
  which routinely happens once block memberships are nearly hard late in a VEM run.
  As a consequence, `nloptr` is no longer a dependency of the package
  (`src/nlopt_wrapper.*`, `src/packing.*` are removed). This is a wash on raw wall-clock
  time in our benchmarks (Newton's per-iteration cost grows with the number of
  covariates, offsetting the much lower iteration count for models with many
  covariates), but removes an external dependency and gives more predictable, reliable
  convergence.
- remove src/utils.h: none of its helpers were actually used anywhere
- minor documentation and dead-code cleanup (see git history for details)

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

