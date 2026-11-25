# Base internal class for adjusting a binary Stochastic Block Model in the context of missSBM.

It is not designed to be called directly by the user; see the concrete
variants
[`SimpleSBM_fit_noCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.md),
[`SimpleSBM_fit_withCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_withCov.md)
and
[`SimpleSBM_fit_MNAR`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_MNAR.md).

## Super classes

[`sbm::SBM`](https://grosssbm.github.io/sbm/reference/SBM.html) -\>
[`sbm::SimpleSBM`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html)
-\> `SimpleSBM_fit`

## Active bindings

- `type`:

  the type of SBM (distribution of edges values, network type, presence
  of covariates)

- `nbObservedDyads`:

  the number of observed dyads

- `penalty`:

  double, value of the penalty term in ICL: the connectivity and the
  effects of the covariates are estimated on the observed dyads, which
  are counted instead of all the possible dyads (as sbm does for a
  network with missing values)

- `entropy`:

  double, value of the entropy due to the clustering distribution

- `vExpec`:

  double: variational expectation of the complete log-likelihood of the
  observed dyads and of the blocks

- `loglik`:

  double: approximation of the log-likelihood (variational lower bound)
  reached

- `ICL`:

  double: the ICL (Integrated Classification Likelihood) of the SBM
  alone: -2 times the variational expectation of the complete
  log-likelihood of the observed dyads at the MAP memberships (each node
  in its most probable block), plus the penalty, as for
  [missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)

## Methods

### Public methods

- [`SimpleSBM_fit$new()`](#method-SimpleSBM_fit-initialize)

- [`SimpleSBM_fit$doVEM()`](#method-SimpleSBM_fit-doVEM)

- [`SimpleSBM_fit$reorder()`](#method-SimpleSBM_fit-reorder)

- [`SimpleSBM_fit$get_state()`](#method-SimpleSBM_fit-get_state)

- [`SimpleSBM_fit$set_state()`](#method-SimpleSBM_fit-set_state)

- [`SimpleSBM_fit$polish_log_tau()`](#method-SimpleSBM_fit-polish_log_tau)

- [`SimpleSBM_fit$clone()`](#method-SimpleSBM_fit-clone)

Inherited methods

- [`sbm::SBM$print()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-print)
- [`sbm::SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)
- [`sbm::SimpleSBM$plot()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-plot)
- [`sbm::SimpleSBM$predict()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-predict)
- [`sbm::SimpleSBM$rEdges()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rEdges)
- [`sbm::SimpleSBM$rMemberships()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rMemberships)
- [`sbm::SimpleSBM$show()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-show)

------------------------------------------------------------------------

### `SimpleSBM_fit$new()`

constructor for simpleSBM_fit for missSBM purpose

#### Usage

    SimpleSBM_fit$new(networkData, clusterInit, covarList = list())

#### Arguments

- `networkData`:

  a structure to store network under missing data condition: either a
  matrix possibly with NA, or a missSBM:::partlyObservedNetwork

- `clusterInit`:

  Initial clustering: a vector with size `ncol(adjacencyMatrix)`,
  providing a user-defined clustering with `nbBlocks` levels.

- `covarList`:

  An optional list with M entries (the M covariates).

------------------------------------------------------------------------

### `SimpleSBM_fit$doVEM()`

method to perform estimation via variational EM

#### Usage

    SimpleSBM_fit$doVEM(
      threshold = 0.01,
      maxIter = 100,
      fixPointIter = 3,
      trace = FALSE
    )

#### Arguments

- `threshold`:

  stop when an optimization step changes the objective function by less
  than threshold. Default is 1e-4.

- `maxIter`:

  V-EM algorithm stops when the number of iteration exceeds maxIter.
  Default is 10

- `fixPointIter`:

  number of fix-point iterations in the Variational E step. Default is
  5.

- `trace`:

  logical for verbosity. Default is `FALSE`.

------------------------------------------------------------------------

### `SimpleSBM_fit$reorder()`

permute group labels by order of decreasing probability

#### Usage

    SimpleSBM_fit$reorder()

#### Returns

invisibly, the permutation `o` applied (new block `k` is former block
`o[k]`), so that the parameters of the sampling indexed by the blocks
can follow

------------------------------------------------------------------------

### `SimpleSBM_fit$get_state()`

a lightweight snapshot of the mutable VEM state (as opposed to
`clone()`, which duplicates the whole object, including the – possibly
large – network data)

#### Usage

    SimpleSBM_fit$get_state()

------------------------------------------------------------------------

### `SimpleSBM_fit$set_state()`

restore a state previously returned by `get_state()`

#### Usage

    SimpleSBM_fit$set_state(state)

#### Arguments

- `state`:

  a state, as returned by `get_state()`

------------------------------------------------------------------------

### `SimpleSBM_fit$polish_log_tau()`

for each node, the log-likelihood of the observed dyads it would
contribute to each block if hard-assigned there (parameters held fixed),
plus the term of the sampling: the gains of the node swaps of
`missSBM_fit$polish()`

#### Usage

    SimpleSBM_fit$polish_log_tau(log_lambda = 0, ...)

#### Arguments

- `log_lambda`:

  term of the sampling depending on the blocks, added as is (0 for the
  designs whose sampling does not depend on them)

- `...`:

  additional arguments, not used

#### Returns

an N x Q matrix

------------------------------------------------------------------------

### `SimpleSBM_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
