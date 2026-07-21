# Class for fitting a block-node sampling

Class for fitting a block-node sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
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

- [`blockNodeSampling_fit$clone()`](#method-blockNodeSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingNodes_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-show)
- [`networkSamplingNodes_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-update_imputation)

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

a method to update the estimation of the parameters. By default, nothing
to do (corresponds to MAR sampling)

#### Usage

    blockNodeSampling_fit$update_parameters(imputedNet, Z)

#### Arguments

- `imputedNet`:

  an adjacency matrix where missing values have been imputed

- `Z`:

  indicator of blocks

------------------------------------------------------------------------

### `blockNodeSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    blockNodeSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
