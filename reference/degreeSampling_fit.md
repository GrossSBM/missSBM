# Class for fitting a degree sampling

Class for fitting a degree sampling

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampling_fit`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.md)
-\>
[`networkSamplingNodes_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.md)
-\> `degreeSampling_fit`

## Active bindings

- `vExpec`:

  lower bound of the variational expectation of the log-likelihood of
  the sampling, by the bound of Jaakkola and Jordan: with x_i = a + b
  D_i, \\\log \sigma(x) \geq \log \sigma(\zeta) + (x - \zeta)/2 +
  h(\zeta)(x^2 - \zeta^2)\\ and \\\log(1 - \sigma(x)) = \log \sigma(x) -
  x\\, in expectation over the imputations

## Methods

### Public methods

- [`degreeSampling_fit$new()`](#method-degreeSampling_fit-initialize)

- [`degreeSampling_fit$update_parameters()`](#method-degreeSampling_fit-update_parameters)

- [`degreeSampling_fit$update_imputation()`](#method-degreeSampling_fit-update_imputation)

- [`degreeSampling_fit$clone()`](#method-degreeSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling_fit$reorder()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-reorder)
- [`networkSampling_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-show)
- [`networkSampling_fit$vExpec_at()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-vExpec_at)

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

update the parameters given the imputations: the expected degrees, the
parameters zeta of the bound of Jaakkola and Jordan, then (a, b)

#### Usage

    degreeSampling_fit$update_parameters(nu, ...)

#### Arguments

- `nu`:

  the imputations of the missing dyads (sparse matrix with the pattern
  of the missing dyads)

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
