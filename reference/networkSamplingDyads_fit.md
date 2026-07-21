# Virtual class used to define a family of networkSamplingDyads_fit

Virtual class used to define a family of networkSamplingDyads_fit

## Super class

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\> `networkSamplingDyads_fit`

## Active bindings

- `penalty`:

  double, value of the penalty term in ICL

- `log_lambda`:

  double, term for adjusting the imputation step which depends on the
  type of sampling

## Methods

### Public methods

- [`networkSamplingDyads_fit$new()`](#method-networkSamplingDyads_fit-initialize)

- [`networkSamplingDyads_fit$show()`](#method-networkSamplingDyads_fit-show)

- [`networkSamplingDyads_fit$update_parameters()`](#method-networkSamplingDyads_fit-update_parameters)

- [`networkSamplingDyads_fit$update_imputation()`](#method-networkSamplingDyads_fit-update_imputation)

- [`networkSamplingDyads_fit$clone()`](#method-networkSamplingDyads_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)

------------------------------------------------------------------------

### `networkSamplingDyads_fit$new()`

constructor for networkSampling_fit

#### Usage

    networkSamplingDyads_fit$new(partlyObservedNetwork, name)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `name`:

  a character for the name of sampling to fit on the
  partlyObservedNetwork

------------------------------------------------------------------------

### `networkSamplingDyads_fit$show()`

show method

#### Usage

    networkSamplingDyads_fit$show()

------------------------------------------------------------------------

### `networkSamplingDyads_fit$update_parameters()`

a method to update the estimation of the parameters. By default, nothing
to do (corresponds to MAR sampling)

#### Usage

    networkSamplingDyads_fit$update_parameters(...)

#### Arguments

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `networkSamplingDyads_fit$update_imputation()`

a method to update the imputation of the missing entries.

#### Usage

    networkSamplingDyads_fit$update_imputation(nu)

#### Arguments

- `nu`:

  the matrix of (uncorrected) imputation for missing entries

------------------------------------------------------------------------

### `networkSamplingDyads_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    networkSamplingDyads_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
