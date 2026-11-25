# Class for fitting a block-node sampling

Class for fitting a block-node sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampling_fit`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.md)
-\>
[`networkSamplingNodes_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.md)
-\> `blockNodeSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

- `log_lambda`:

  double, term for adjusting the imputation step which depends on the
  type of sampling

## Methods

### Public methods

- [`blockNodeSampling_fit$new()`](#method-blockNodeSampling_fit-initialize)

- [`blockNodeSampling_fit$update_parameters()`](#method-blockNodeSampling_fit-update_parameters)

- [`blockNodeSampling_fit$reorder()`](#method-blockNodeSampling_fit-reorder)

- [`blockNodeSampling_fit$vExpec_at()`](#method-blockNodeSampling_fit-vExpec_at)

- [`blockNodeSampling_fit$clone()`](#method-blockNodeSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-show)
- [`networkSampling_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-update_imputation)

------------------------------------------------------------------------

### `blockNodeSampling_fit$new()`

constructor

#### Usage

    blockNodeSampling_fit$new(partlyObservedNetwork, blockInit)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `blockInit`:

  n x Q matrix of initial block indicators

------------------------------------------------------------------------

### `blockNodeSampling_fit$update_parameters()`

update the parameters: the observation rate of the nodes of each block

#### Usage

    blockNodeSampling_fit$update_parameters(Z, ...)

#### Arguments

- `Z`:

  probabilities of block memberships

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `blockNodeSampling_fit$reorder()`

permute the parameters indexed by the blocks, after the blocks of the
SBM were permuted

#### Usage

    blockNodeSampling_fit$reorder(o)

#### Arguments

- `o`:

  the permutation applied to the blocks (new block k is former block
  `o[k]`)

------------------------------------------------------------------------

### `blockNodeSampling_fit$vExpec_at()`

the variational expectation of the log-likelihood of the sampling at
given memberships (used by the ICL, at the MAP memberships)

#### Usage

    blockNodeSampling_fit$vExpec_at(Z)

#### Arguments

- `Z`:

  an n x Q matrix of memberships

------------------------------------------------------------------------

### `blockNodeSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    blockNodeSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
