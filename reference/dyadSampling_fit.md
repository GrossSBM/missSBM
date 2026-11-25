# Class for fitting a dyad sampling

Class for fitting a dyad sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampling_fit`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.md)
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
- [`networkSampling_fit$reorder()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-reorder)
- [`networkSampling_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-show)
- [`networkSampling_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-update_imputation)
- [`networkSampling_fit$update_parameters()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-update_parameters)
- [`networkSampling_fit$vExpec_at()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-vExpec_at)

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
