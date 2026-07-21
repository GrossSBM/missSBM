# Internal class for a binary SBM fit under MNAR sampling designs (double-standard, block-node, block-dyad).

It is not designed to be called directly by the user.

## Super classes

[`sbm::SBM`](https://grosssbm.github.io/sbm/reference/SBM.html) -\>
[`sbm::SimpleSBM`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html)
-\>
[`SimpleSBM_fit`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.md)
-\>
[`SimpleSBM_fit_noCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.md)
-\> `SimpleSBM_MNAR_noCov`

## Active bindings

- `vExpec`:

  double: variational approximation of the expectation complete
  log-likelihood

## Methods

### Public methods

- [`SimpleSBM_MNAR_noCov$new()`](#method-SimpleSBM_MNAR_noCov-initialize)

- [`SimpleSBM_MNAR_noCov$update_parameters()`](#method-SimpleSBM_MNAR_noCov-update_parameters)

- [`SimpleSBM_MNAR_noCov$update_blocks()`](#method-SimpleSBM_MNAR_noCov-update_blocks)

- [`SimpleSBM_MNAR_noCov$polish_log_tau()`](#method-SimpleSBM_MNAR_noCov-polish_log_tau)

- [`SimpleSBM_MNAR_noCov$clone()`](#method-SimpleSBM_MNAR_noCov-clone)

Inherited methods

- [`sbm::SBM$print()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-print)
- [`sbm::SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)
- [`sbm::SimpleSBM$plot()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-plot)
- [`sbm::SimpleSBM$predict()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-predict)
- [`sbm::SimpleSBM$rEdges()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rEdges)
- [`sbm::SimpleSBM$rMemberships()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rMemberships)
- [`sbm::SimpleSBM$show()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-show)
- [`SimpleSBM_fit$doVEM()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-doVEM)
- [`SimpleSBM_fit$get_state()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-get_state)
- [`SimpleSBM_fit$reorder()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-reorder)
- [`SimpleSBM_fit$set_state()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-set_state)

------------------------------------------------------------------------

### `SimpleSBM_MNAR_noCov$new()`

constructor for simpleSBM_fit for missSBM purpose

#### Usage

    SimpleSBM_MNAR_noCov$new(networkData, clusterInit)

#### Arguments

- `networkData`:

  a structure to store network under missing data condition: either a
  matrix possibly with NA, or a missSBM:::partlyObservedNetwork

- `clusterInit`:

  Initial clustering: a vector with size `ncol(adjacencyMatrix)`,
  providing a user-defined clustering with `nbBlocks` levels.

------------------------------------------------------------------------

### `SimpleSBM_MNAR_noCov$update_parameters()`

update parameters estimation (M-step)

#### Usage

    SimpleSBM_MNAR_noCov$update_parameters(nu = NULL)

#### Arguments

- `nu`:

  currently imputed values

------------------------------------------------------------------------

### `SimpleSBM_MNAR_noCov$update_blocks()`

update variational estimation of blocks (VE-step)

#### Usage

    SimpleSBM_MNAR_noCov$update_blocks(log_lambda = 0)

#### Arguments

- `log_lambda`:

  additional term sampling dependent used to de-bias estimation of tau

------------------------------------------------------------------------

### `SimpleSBM_MNAR_noCov$polish_log_tau()`

for each node, the complete-data log-likelihood it would contribute to
each class if hard-assigned there (theta/pi held fixed), used to decide
node-swap moves in `missSBM_fit$polish()`. Unlike `update_blocks()`'s
tau update (which combines the observed- and imputed-part E-steps
as-is), this subtracts one copy of `log(pi)`: both E-step calls add it,
so summing them double-counts it relative to `vExpec`'s convention –
verified numerically to matter (see git history).

#### Usage

    SimpleSBM_MNAR_noCov$polish_log_tau(log_lambda = 0)

#### Arguments

- `log_lambda`:

  additional sampling-design-dependent term, added as-is

#### Returns

an N x Q matrix

------------------------------------------------------------------------

### `SimpleSBM_MNAR_noCov$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_MNAR_noCov$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
