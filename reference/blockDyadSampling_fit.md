# Class for fitting a block-dyad sampling

Class for fitting a block-dyad sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSamplingDyads_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.md)
-\> `blockDyadSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

- `log_lambda`:

  matrix, term for adjusting the imputation step which depends on the
  type of sampling

## Methods

### Public methods

- [`blockDyadSampling_fit$new()`](#method-blockDyadSampling_fit-initialize)

- [`blockDyadSampling_fit$update_parameters()`](#method-blockDyadSampling_fit-update_parameters)

- [`blockDyadSampling_fit$clone()`](#method-blockDyadSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingDyads_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-show)
- [`networkSamplingDyads_fit$update_imputation()`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.html#method-update_imputation)

------------------------------------------------------------------------

### `blockDyadSampling_fit$new()`

constructor

#### Usage

    blockDyadSampling_fit$new(partlyObservedNetwork, blockInit)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `blockInit`:

  n x Q matrix of initial block indicators

------------------------------------------------------------------------

### `blockDyadSampling_fit$update_parameters()`

a method to update the estimation of the parameters. By default, nothing
to do (corresponds to MAR sampling)

#### Usage

    blockDyadSampling_fit$update_parameters(nu, Z)

#### Arguments

- `nu`:

  the matrix of (uncorrected) imputation for missing entries

- `Z`:

  probabilities of block memberships

------------------------------------------------------------------------

### `blockDyadSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    blockDyadSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
