# Class for defining a block dyad sampler

Class for defining a block dyad sampler

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampler`](https://grosssbm.github.io/missSBM/reference/networkSampler.md)
-\>
[`dyadSampler`](https://grosssbm.github.io/missSBM/reference/dyadSampler.md)
-\> `blockDyadSampler`

## Active bindings

- `df`:

  the number of parameters of this sampling

## Methods

### Public methods

- [`blockDyadSampler$new()`](#method-blockDyadSampler-initialize)

- [`blockDyadSampler$clone()`](#method-blockDyadSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)
- [`networkSampler$rSamplingMatrix()`](https://grosssbm.github.io/missSBM/reference/networkSampler.html#method-rSamplingMatrix)

------------------------------------------------------------------------

### `blockDyadSampler$new()`

constructor for networkSampling

#### Usage

    blockDyadSampler$new(
      parameters = NA,
      nbNodes = NA,
      directed = FALSE,
      clusters = NA
    )

#### Arguments

- `parameters`:

  the vector of parameters associated to the sampling at play

- `nbNodes`:

  number of nodes in the network

- `directed`:

  logical, directed network of not

- `clusters`:

  a vector of class memberships

------------------------------------------------------------------------

### `blockDyadSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    blockDyadSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
