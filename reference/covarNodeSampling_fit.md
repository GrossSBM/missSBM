# Class for fitting a node-centered sampling with covariate

Class for fitting a node-centered sampling with covariate

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSamplingNodes_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.md)
-\> `covarNodeSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

## Methods

### Public methods

- [`covarNodeSampling_fit$new()`](#method-covarNodeSampling_fit-initialize)

- [`covarNodeSampling_fit$clone()`](#method-covarNodeSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingNodes_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-show)
- [`networkSamplingNodes_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-update_imputation)
- [`networkSamplingNodes_fit$update_parameters()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-update_parameters)

------------------------------------------------------------------------

### `covarNodeSampling_fit$new()`

constructor

#### Usage

    covarNodeSampling_fit$new(partlyObservedNetwork, ...)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `...`:

  used for compatibility

------------------------------------------------------------------------

### `covarNodeSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    covarNodeSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
