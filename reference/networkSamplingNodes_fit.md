# Virtual class used to define a family of networkSamplingNodes_fit

Virtual class used to define a family of networkSamplingNodes_fit

## Super class

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\> `networkSamplingNodes_fit`

## Active bindings

- `penalty`:

  double, value of the penalty term in ICL

- `log_lambda`:

  double, term for adjusting the imputation step which depends on the
  type of sampling

## Methods

### Public methods

- [`networkSamplingNodes_fit$new()`](#method-networkSamplingNodes_fit-initialize)

- [`networkSamplingNodes_fit$show()`](#method-networkSamplingNodes_fit-show)

- [`networkSamplingNodes_fit$update_parameters()`](#method-networkSamplingNodes_fit-update_parameters)

- [`networkSamplingNodes_fit$update_imputation()`](#method-networkSamplingNodes_fit-update_imputation)

- [`networkSamplingNodes_fit$clone()`](#method-networkSamplingNodes_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)

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

### `networkSamplingNodes_fit$show()`

show method

#### Usage

    networkSamplingNodes_fit$show()

------------------------------------------------------------------------

### `networkSamplingNodes_fit$update_parameters()`

a method to update the estimation of the parameters. By default, nothing
to do (corresponds to MAR sampling)

#### Usage

    networkSamplingNodes_fit$update_parameters(...)

#### Arguments

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `networkSamplingNodes_fit$update_imputation()`

a method to update the imputation of the missing entries.

#### Usage

    networkSamplingNodes_fit$update_imputation(nu)

#### Arguments

- `nu`:

  the matrix of (uncorrected) imputation for missing entries

------------------------------------------------------------------------

### `networkSamplingNodes_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    networkSamplingNodes_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
