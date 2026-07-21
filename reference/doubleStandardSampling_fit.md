# Class for fitting a double-standard sampling

Class for fitting a double-standard sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSamplingDyads_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.md)
-\> `doubleStandardSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

## Methods

### Public methods

- [`doubleStandardSampling_fit$new()`](#method-doubleStandardSampling_fit-initialize)

- [`doubleStandardSampling_fit$update_parameters()`](#method-doubleStandardSampling_fit-update_parameters)

- [`doubleStandardSampling_fit$update_imputation()`](#method-doubleStandardSampling_fit-update_imputation)

- [`doubleStandardSampling_fit$clone()`](#method-doubleStandardSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingDyads_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-show)

------------------------------------------------------------------------

### `doubleStandardSampling_fit$new()`

constructor

#### Usage

    doubleStandardSampling_fit$new(partlyObservedNetwork, ...)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `...`:

  used for compatibility

------------------------------------------------------------------------

### `doubleStandardSampling_fit$update_parameters()`

a method to update the estimation of the parameters. By default, nothing
to do (corresponds to MAR sampling)

#### Usage

    doubleStandardSampling_fit$update_parameters(nu, ...)

#### Arguments

- `nu`:

  an adjacency matrix with imputed values (only)

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `doubleStandardSampling_fit$update_imputation()`

a method to update the imputation of the missing entries.

#### Usage

    doubleStandardSampling_fit$update_imputation(nu)

#### Arguments

- `nu`:

  the matrix of (uncorrected) imputation for missing entries

------------------------------------------------------------------------

### `doubleStandardSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    doubleStandardSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
