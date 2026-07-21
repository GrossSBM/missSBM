# Class for defining a simple dyad sampler

Class for defining a simple dyad sampler

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampler`](https://grosssbm.github.io/missSBM/reference/networkSampler.md)
-\>
[`dyadSampler`](https://grosssbm.github.io/missSBM/reference/dyadSampler.md)
-\> `simpleDyadSampler`

## Methods

### Public methods

- [`simpleDyadSampler$new()`](#method-simpleDyadSampler-initialize)

- [`simpleDyadSampler$clone()`](#method-simpleDyadSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)
- [`networkSampler$rSamplingMatrix()`](https://grosssbm.github.io/missSBM/reference/networkSampler.html#method-rSamplingMatrix)

------------------------------------------------------------------------

### `simpleDyadSampler$new()`

constructor for networkSampling

#### Usage

    simpleDyadSampler$new(
      parameters = NA,
      nbNodes = NA,
      directed = FALSE,
      covarArray = NULL,
      intercept = 0
    )

#### Arguments

- `parameters`:

  the vector of parameters associated to the sampling at play

- `nbNodes`:

  number of nodes in the network

- `directed`:

  logical, directed network of not

- `covarArray`:

  an array of covariates used

- `intercept`:

  double, intercept term used to compute the probability of sampling in
  the presence of covariates. Default 0.

------------------------------------------------------------------------

### `simpleDyadSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    simpleDyadSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
