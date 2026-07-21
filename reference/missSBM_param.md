# Control of a missSBM fit

Helper to define the list of parameters that controls
[`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md).
All arguments have defaults.

## Usage

``` r
missSBM_param(
  threshold = 0.01,
  maxIter = 50,
  fixPointIter = 3,
  imputation = c("median", "average", "zero"),
  similarity = l1_similarity,
  useCov = TRUE,
  clusterInit = NULL,
  polish = TRUE,
  iterates = 1,
  maxMergeCandidates = 30,
  stopOnDegenerate = TRUE,
  maxConsecutiveDegenerate = 2,
  warmChain = FALSE,
  trace = TRUE
)
```

## Arguments

- threshold:

  V-EM algorithm stops when an optimization step changes the objective
  function or the parameters by less than threshold. Default is 1e-2.

- maxIter:

  V-EM algorithm stops when the number of iteration exceeds maxIter.
  Default is 50.

- fixPointIter:

  number of fix-point iterations in the V-E step. Default is 3.

- imputation:

  character, the imputation strategy used to build the initial
  clustering (see
  [`partlyObservedNetwork`](https://grosssbm.github.io/missSBM/reference/partlyObservedNetwork.md)'s
  `imputation()`). Either "median" (default), "average" or "zero".

- similarity:

  an R x R -\> R function to compute similarities between node
  covariates. Default is
  [`l1_similarity`](https://grosssbm.github.io/missSBM/reference/l1_similarity.md),
  that is, -abs(x - y). Only relevant when the covariates are
  node-centered (i.e. `covariates` is a list of size-N vectors).

- useCov:

  logical. If `covariates` is not empty, should they be used for the SBM
  inference (or just for the sampling)? Default is TRUE.

- clusterInit:

  initial clustering: `NULL` (default) for a cold spectral clustering,
  or a list with `length(vBlocks)` vectors, each with size
  `ncol(adjacencyMatrix)`, providing a user-defined clustering.

- polish:

  logical, should each model be node-swap-polished (see
  [missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
  `polish()`) after the initial VEM fit, fixing individually
  misclassified nodes at that model's own number of blocks? Cheap
  relative to exploration (`iterates`), since it does not search across
  numbers of blocks. Default is TRUE.

- iterates:

  integer, the number of forward/backward exploration passes searching
  for a better number of blocks by splitting/merging clusters (unlike
  `polish`, which only refines each model at its own number of blocks):
  more expensive than `polish`. `0` disables exploration entirely (only
  the initial VEM fit, plus `polish` if enabled, is returned). Default
  is 1.

- maxMergeCandidates:

  integer, caps the number of cluster-pair merge candidates tried during
  backward exploration (quadratic in the number of blocks otherwise).
  Beyond this cap, only the pairs with the most similar fitted
  connectivity profiles are tried, since merging two blocks with very
  different connectivity is rarely competitive anyway. Default is 30.

- stopOnDegenerate:

  logical. A requested number of blocks can be higher than what the
  network actually supports: VEM then collapses one or more classes (see
  [missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
  `repair()`), and even a fair recovery attempt (`repair()`, then a full
  exploration pass letting `explore_forward()` try to split a healthy
  neighbor into that range) can still collapse right back down. When
  this persists over `maxConsecutiveDegenerate` consecutive (increasing)
  values of `vBlocks` after a pass, forward (split) exploration stops
  growing further into that range on subsequent passes (only relevant
  when `iterates > 1`); the affected models stay as fitted, flagged via
  `$degenerate`, and a warning is issued. Default is TRUE.

- maxConsecutiveDegenerate:

  integer, the run length that triggers `stopOnDegenerate`. Default is
  2.

- warmChain:

  logical (**work in progress**). If `TRUE`, each model is initialized
  by splitting (see
  [missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
  [`split()`](https://rdrr.io/r/base/split.html)) the already-converged,
  smaller neighbor in `vBlocks` instead of an independent cold spectral
  clustering (see
  [missSBM_collection](https://grosssbm.github.io/missSBM/reference/missSBM_collection.md)'s
  `estimate_chain()`) – meant to reduce VEM component collapse at higher
  numbers of blocks. Sequential in the number of blocks, unlike the
  default (`FALSE`), which fits every model in parallel: can be slower
  in wall-clock time when several workers are available. Default is
  FALSE.

- trace:

  logical for verbosity. Default is TRUE.

## Value

a list of parameters configuring the fit, with class `missSBM_param`.

## See also

[`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)

## Examples

``` r
my_control <- missSBM_param(iterates = 2, polish = FALSE)
my_control$iterates
#> [1] 2
```
