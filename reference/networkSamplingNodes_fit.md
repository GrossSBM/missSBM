# Virtual class of the fits of the sampling designs on the nodes

Virtual class of the fits of the sampling designs on the nodes

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampling_fit`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.md)
-\> `networkSamplingNodes_fit`

## Methods

### Public methods

- [`networkSamplingNodes_fit$new()`](#method-networkSamplingNodes_fit-initialize)

- [`networkSamplingNodes_fit$clone()`](#method-networkSamplingNodes_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling_fit$reorder()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-reorder)
- [`networkSampling_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-show)
- [`networkSampling_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-update_imputation)
- [`networkSampling_fit$update_parameters()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-update_parameters)
- [`networkSampling_fit$vExpec_at()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-vExpec_at)

------------------------------------------------------------------------

### `networkSamplingNodes_fit$new()`

constructor

#### Usage

    networkSamplingNodes_fit$new(partlyObservedNetwork, name)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `name`:

  a character for the name of sampling to fit on the
  partlyObservedNetwork

------------------------------------------------------------------------

### `networkSamplingNodes_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    networkSamplingNodes_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
