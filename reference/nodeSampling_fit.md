# Class for fitting a node sampling

Class for fitting a node sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSamplingNodes_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.md)
-\> `nodeSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

## Methods

### Public methods

- [`nodeSampling_fit$new()`](#method-nodeSampling_fit-initialize)

- [`nodeSampling_fit$clone()`](#method-nodeSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingNodes_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-show)
- [`networkSamplingNodes_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-update_imputation)
- [`networkSamplingNodes_fit$update_parameters()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-update_parameters)

------------------------------------------------------------------------

### `nodeSampling_fit$new()`

constructor

#### Usage

    nodeSampling_fit$new(partlyObservedNetwork, ...)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `...`:

  used for compatibility

------------------------------------------------------------------------

### `nodeSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    nodeSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
