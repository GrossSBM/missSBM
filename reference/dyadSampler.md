# Virtual class for all dyad-centered samplers

Virtual class for all dyad-centered samplers

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampler`](https://grosssbm.github.io/missSBM/reference/networkSampler.md)
-\> `dyadSampler`

## Methods

### Public methods

- [`dyadSampler$new()`](#method-dyadSampler-initialize)

- [`dyadSampler$clone()`](#method-dyadSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)
- [`networkSampler$rSamplingMatrix()`](https://grosssbm.github.io/missSBM/reference/networkSampler.html#method-rSamplingMatrix)

------------------------------------------------------------------------

### `dyadSampler$new()`

constructor for networkSampling

#### Usage

    dyadSampler$new(type = NA, parameters = NA, nbNodes = NA, directed = FALSE)

#### Arguments

- `type`:

  character for the type of sampling. must be in ("dyad", "covar-dyad",
  "node", "covar-node", "block-node", "block-dyad", "double-standard",
  "degree")

- `parameters`:

  the vector of parameters associated to the sampling at play

- `nbNodes`:

  number of nodes in the network

- `directed`:

  logical, directed network of not

------------------------------------------------------------------------

### `dyadSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    dyadSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
