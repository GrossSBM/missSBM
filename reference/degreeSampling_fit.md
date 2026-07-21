# Class for fitting a degree sampling

Class for fitting a degree sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSamplingNodes_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.md)
-\> `degreeSampling_fit`

## Active bindings

- `vExpec`:

  variational expectation of the sampling

## Methods

### Public methods

- [`degreeSampling_fit$new()`](#method-degreeSampling_fit-initialize)

- [`degreeSampling_fit$update_parameters()`](#method-degreeSampling_fit-update_parameters)

- [`degreeSampling_fit$update_imputation()`](#method-degreeSampling_fit-update_imputation)

- [`degreeSampling_fit$clone()`](#method-degreeSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSamplingNodes_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.html#method-show)

------------------------------------------------------------------------

### `degreeSampling_fit$new()`

constructor

#### Usage

    degreeSampling_fit$new(partlyObservedNetwork, blockInit, connectInit)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `blockInit`:

  n x Q matrix of initial block indicators

- `connectInit`:

  Q x Q matrix of initial block probabilities of connection

------------------------------------------------------------------------

### `degreeSampling_fit$update_parameters()`

a method to update the estimation of the parameters. By default, nothing
to do (corresponds to MAR sampling)

#### Usage

    degreeSampling_fit$update_parameters(imputedNet, ...)

#### Arguments

- `imputedNet`:

  an adjacency matrix where missing values have been imputed

- `...`:

  used for compatibility

------------------------------------------------------------------------

### `degreeSampling_fit$update_imputation()`

a method to update the imputation of the missing entries.

#### Usage

    degreeSampling_fit$update_imputation(PI, ...)

#### Arguments

- `PI`:

  the matrix of inter/intra class probability of connection

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `degreeSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    degreeSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
