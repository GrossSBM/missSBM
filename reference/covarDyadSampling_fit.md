# Class for fitting a dyad sampling with covariates

Class for fitting a dyad sampling with covariates

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSamplingDyads_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.md)
-\> `covarDyadSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

## Methods

### Public methods

- [`covarDyadSampling_fit$new()`](#method-covarDyadSampling_fit-initialize)

- [`covarDyadSampling_fit$clone()`](#method-covarDyadSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingDyads_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-show)
- [`networkSamplingDyads_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-update_imputation)
- [`networkSamplingDyads_fit$update_parameters()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-update_parameters)

------------------------------------------------------------------------

### `covarDyadSampling_fit$new()`

constructor

#### Usage

    covarDyadSampling_fit$new(partialNet, ...)

#### Arguments

- `partialNet`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `...`:

  used for compatibility

------------------------------------------------------------------------

### `covarDyadSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    covarDyadSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
