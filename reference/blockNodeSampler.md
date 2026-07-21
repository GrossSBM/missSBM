# Class for defining a block node sampler

Class for defining a block node sampler

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampler`](https://grosssbm.github.io/missSBM/reference/networkSampler.md)
-\>
[`nodeSampler`](https://grosssbm.github.io/missSBM/reference/nodeSampler.md)
-\> `blockNodeSampler`

## Methods

### Public methods

- [`blockNodeSampler$new()`](#method-blockNodeSampler-initialize)

- [`blockNodeSampler$clone()`](#method-blockNodeSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)
- [`networkSampler$rSamplingMatrix()`](https://grosssbm.github.io/missSBM/reference/networkSampler.html#method-rSamplingMatrix)

------------------------------------------------------------------------

### `blockNodeSampler$new()`

constructor for networkSampling

#### Usage

    blockNodeSampler$new(
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

### `blockNodeSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    blockNodeSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
