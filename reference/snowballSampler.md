# Class for defining a snowball sampler

Class for defining a snowball sampler

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampler`](https://grosssbm.github.io/missSBM/reference/networkSampler.md)
-\>
[`nodeSampler`](https://grosssbm.github.io/missSBM/reference/nodeSampler.md)
-\> `snowballSampler`

## Methods

### Public methods

- [`snowballSampler$new()`](#method-snowballSampler-initialize)

- [`snowballSampler$clone()`](#method-snowballSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)
- [`networkSampler$rSamplingMatrix()`](https://grosssbm.github.io/missSBM/reference/networkSampler.html#method-rSamplingMatrix)

------------------------------------------------------------------------

### `snowballSampler$new()`

constructor for networkSampling

#### Usage

    snowballSampler$new(parameters = NA, adjacencyMatrix = NA, directed = FALSE)

#### Arguments

- `parameters`:

  the vector of parameters associated to the sampling at play

- `adjacencyMatrix`:

  the adjacency matrix of the network

- `directed`:

  logical, directed network of not

------------------------------------------------------------------------

### `snowballSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    snowballSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
