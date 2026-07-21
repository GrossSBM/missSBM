# Class for fitting a dyad sampling

Class for fitting a dyad sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSamplingDyads_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.md)
-\> `dyadSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

## Methods

### Public methods

- [`dyadSampling_fit$new()`](#method-dyadSampling_fit-initialize)

- [`dyadSampling_fit$clone()`](#method-dyadSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingDyads_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-show)
- [`networkSamplingDyads_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-update_imputation)
- [`networkSamplingDyads_fit$update_parameters()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-update_parameters)

------------------------------------------------------------------------

### `dyadSampling_fit$new()`

constructor

#### Usage

    dyadSampling_fit$new(partlyObservedNetwork, ...)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `...`:

  used for compatibility

------------------------------------------------------------------------

### `dyadSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    dyadSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
