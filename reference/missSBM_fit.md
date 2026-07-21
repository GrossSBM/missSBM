# An R6 class to represent an SBM fit with missing data

The function
[`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
fits a collection of SBM for varying number of block. Each fitted SBM is
an instance of an R6 object with class `missSBM_fit`, described here.

Fields are accessed via active binding and cannot be changed by the
user.

This class comes with a set of R6 methods, some of them being useful for
the user and exported as S3 methods. See the documentation for
[`show()`](https://rdrr.io/r/methods/show.html),
[`print()`](https://rdrr.io/r/base/print.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Active bindings

- `fittedSBM`:

  the fitted SBM with class
  [`SimpleSBM_fit_noCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.md),
  [`SimpleSBM_fit_withCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_withCov.md)
  or
  [`SimpleSBM_fit_MNAR`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_MNAR.md)
  inheriting from class
  [`sbm::SimpleSBM_fit`](https://grosssbm.github.io/sbm/reference/SimpleSBM_fit.html)

- `fittedSampling`:

  the fitted sampling, inheriting from class
  [`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
  and corresponding fits

- `imputedNetwork`:

  The network data as a matrix with NAs values imputed with the current
  model

- `monitoring`:

  a list carrying information about the optimization process

- `occupiedBlocks`:

  the number of classes actually occupied by at least one node; can be
  less than `fittedSBM$nbBlocks` for an over-specified fit whose VEM has
  collapsed one or more classes (see `repair()`)

- `degenerate`:

  `TRUE` if `occupiedBlocks < fittedSBM$nbBlocks`

- `entropyImputed`:

  the entropy of the distribution of the imputed dyads

- `entropy`:

  the entropy due to the distribution of the imputed dyads and of the
  clustering

- `vExpec`:

  double: variational expectation of the complete log-likelihood

- `penalty`:

  double, value of the penalty term in ICL

- `loglik`:

  double: approximation of the log-likelihood (variational lower bound)
  reached

- `ICL`:

  double: value of the integrated classification log-likelihood

## Methods

### Public methods

- [`missSBM_fit$new()`](#method-missSBM_fit-initialize)

- [`missSBM_fit$doVEM()`](#method-missSBM_fit-doVEM)

- [`missSBM_fit$split()`](#method-missSBM_fit-split)

- [`missSBM_fit$candidates_split()`](#method-missSBM_fit-candidates_split)

- [`missSBM_fit$merge()`](#method-missSBM_fit-merge)

- [`missSBM_fit$candidates_merge()`](#method-missSBM_fit-candidates_merge)

- [`missSBM_fit$repair()`](#method-missSBM_fit-repair)

- [`missSBM_fit$polish()`](#method-missSBM_fit-polish)

- [`missSBM_fit$show()`](#method-missSBM_fit-show)

- [`missSBM_fit$print()`](#method-missSBM_fit-print)

- [`missSBM_fit$clone()`](#method-missSBM_fit-clone)

------------------------------------------------------------------------

### `missSBM_fit$new()`

constructor for networkSampling

#### Usage

    missSBM_fit$new(partlyObservedNet, netSampling, clusterInit, useCov = TRUE)

#### Arguments

- `partlyObservedNet`:

  An object with class
  [`partlyObservedNetwork`](https://grosssbm.github.io/missSBM/reference/partlyObservedNetwork.md).

- `netSampling`:

  The sampling design for the modelling of missing data: MAR designs
  ("dyad", "node") and MNAR designs ("double-standard", "block-dyad",
  "block-node" ,"degree")

- `clusterInit`:

  Initial clustering: a vector with size `ncol(adjacencyMatrix)`,
  providing a user-defined clustering. The number of blocks is deduced
  from the number of levels in with `clusterInit`.

- `useCov`:

  logical. If covariates are present in partlyObservedNet, should they
  be used for the inference or of the network sampling design, or just
  for the SBM inference? default is TRUE.

------------------------------------------------------------------------

### `missSBM_fit$doVEM()`

a method to perform inference of the current missSBM fit with
variational EM

#### Usage

    missSBM_fit$doVEM(
      control = list(threshold = 0.01, maxIter = 100, fixPointIter = 3, trace = TRUE)
    )

#### Arguments

- `control`:

  a list of VEM control parameters (see
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md))

------------------------------------------------------------------------

### `missSBM_fit$split()`

clone of the current fit after splitting cluster `index` in two, via a
spectral bipartition of the sub-network it induces. Builds but does not
fit the candidate (see `candidates_split()`).

#### Usage

    missSBM_fit$split(index, in_place = FALSE, base_net = NULL)

#### Arguments

- `index`:

  index (integer) of the cluster to split

- `in_place`:

  replace `self`'s own fit (`TRUE`) or return a new object (`FALSE`, the
  default)?

- `base_net`:

  optional precomputed network to bipartition (as built internally at
  the top of this method); lets `candidates_split()` avoid recomputing
  it once per candidate.

#### Returns

a new `missSBM_fit` with one more block, or `NULL` if `index` cannot be
split (its induced sub-network has zero variance)

------------------------------------------------------------------------

### `missSBM_fit$candidates_split()`

generate and cheaply trial-fit candidates obtained by splitting each
splittable cluster in two (see
[`split()`](https://rdrr.io/r/base/split.html)). A cluster is splittable
if it has at least 4 members and non-zero variance in its induced
sub-network.

#### Usage

    missSBM_fit$candidates_split(
      control = list(threshold = 0.01, maxIter = 100, fixPointIter = 3, trace = TRUE),
      trial_niter = 2
    )

#### Arguments

- `control`:

  a list of VEM control parameters (see
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md));
  `maxIter` is overridden by `trial_niter`

- `trial_niter`:

  number of VEM iterations used for the trial fits. Default is 2.

#### Returns

a list of trial-fitted `missSBM_fit` candidates (one per splittable
cluster)

------------------------------------------------------------------------

### `missSBM_fit$merge()`

clone of the current fit after merging clusters `indices[1]` and
`indices[2]` into one. Builds but does not fit the candidate (see
`candidates_merge()`).

#### Usage

    missSBM_fit$merge(indices, in_place = FALSE)

#### Arguments

- `indices`:

  indices (couple of integers) of the clusters to merge

- `in_place`:

  replace `self`'s own fit (`TRUE`) or return a new object (`FALSE`, the
  default)?

#### Returns

a new `missSBM_fit` with one fewer block

------------------------------------------------------------------------

### `missSBM_fit$candidates_merge()`

generate and cheaply trial-fit candidates obtained by merging pairs of
clusters (see [`merge()`](https://rdrr.io/r/base/merge.html)). Beyond
`max_candidates` pairs (quadratic in the number of blocks), only the
most similar-connectivity pairs are tried.

#### Usage

    missSBM_fit$candidates_merge(
      control = list(threshold = 0.01, maxIter = 100, fixPointIter = 3, trace = TRUE),
      max_candidates = 30,
      trial_niter = 2
    )

#### Arguments

- `control`:

  a list of VEM control parameters (see
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md));
  `maxIter` is overridden by `trial_niter`

- `max_candidates`:

  cap on the number of pairs tried. Default is 30.

- `trial_niter`:

  number of VEM iterations used for the trial fits. Default is 2.

#### Returns

a list of trial-fitted `missSBM_fit` candidates

------------------------------------------------------------------------

### `missSBM_fit$repair()`

recovers a degenerate fit (fewer occupied classes than its structural
`nbBlocks`, e.g. after a VEM component collapse) by filling the empty
classes (see `repair_empty_classes()`) and refitting the full VEM.
Mutates `self` in place; a no-op if the fit is not degenerate.

#### Usage

    missSBM_fit$repair(
      control = list(threshold = 0.01, maxIter = 100, fixPointIter = 3, trace = TRUE)
    )

#### Arguments

- `control`:

  a list of VEM control parameters (see
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md))

#### Returns

invisibly, `self`

------------------------------------------------------------------------

### `missSBM_fit$polish()`

discrete node-swap polishing (Kernighan-Lin / greedy-ICL style): after
VEM convergence, tau is near-hard and its fixed point cannot relocate a
single misclassified node (only
[`split()`](https://rdrr.io/r/base/split.html)/[`merge()`](https://rdrr.io/r/base/merge.html)
fix group-level mistakes). Each sweep computes, for every node, the
closed-form complete-data log-likelihood gain of moving it to its best
alternative class (theta/pi held fixed), applies the improving,
non-class-emptying moves, then runs a full VEM to resettle. Stops as
soon as a sweep fails to improve the ICL, mutates `self` in place, and
never leaves the ICL worse than before the call.

#### Usage

    missSBM_fit$polish(
      control = list(threshold = 0.01, maxIter = 100, fixPointIter = 3, trace = TRUE),
      max_sweeps = 10
    )

#### Arguments

- `control`:

  a list of VEM control parameters (see
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md))

- `max_sweeps`:

  maximum number of swap sweeps. Default is 10.

#### Returns

invisibly, `self`

------------------------------------------------------------------------

### `missSBM_fit$show()`

show method for missSBM_fit

#### Usage

    missSBM_fit$show()

------------------------------------------------------------------------

### `missSBM_fit$print()`

User friendly print method

#### Usage

    missSBM_fit$print()

------------------------------------------------------------------------

### `missSBM_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    missSBM_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## Sample 75% of dyads in  French political Blogosphere's network data
adjMatrix <- missSBM::frenchblog2007 %>%
  igraph::as_adjacency_matrix(sparse = FALSE) %>%
  missSBM::observeNetwork(sampling = "dyad", parameters = 0.75)
collection <- estimateMissSBM(adjMatrix, 3:5, sampling = "dyad")
#> 
#> 
#>  Adjusting Variational EM for Stochastic Block Model
#> 
#>  Imputation assumes a 'dyad' network-sampling process
#> 
#>  Initialization of 3 model(s). 
#>  Performing VEM inference
#>      Model with 3 blocks.    Model with 4 blocks.    Model with 5 blocks. Polishing (node-swap)
#> 
#>  Looking for better solutions
#>  Pass 1   Going forward ++                                                                                                     Pass 1   Going backward ++                                                                                                    
my_missSBM_fit <- collection$bestModel
class(my_missSBM_fit)
#> [1] "missSBM_fit" "R6"         
plot(my_missSBM_fit, "imputed")

```
