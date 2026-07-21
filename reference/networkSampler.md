# Definition of R6 Class 'networkSampling_sampler'

This class is use to define a sampling model for a network. Inherits
from 'networkSampling'. Owns a rSampling method which takes an adjacency
matrix as an input and send back an object with class
partlyObservedNetwork.

## See also

[`partlyObservedNetwork`](https://grosssbm.github.io/missSBM/reference/partlyObservedNetwork.md)

## Super class

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\> `networkSampler`

## Active bindings

- `samplingMatrix`:

  a matrix of logical indicating observed entries

## Methods

### Public methods

- [`networkSampler$new()`](#method-networkSampler-initialize)

- [`networkSampler$rSamplingMatrix()`](#method-networkSampler-rSamplingMatrix)

- [`networkSampler$clone()`](#method-networkSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)

------------------------------------------------------------------------

### `networkSampler$new()`

constructor for networkSampling

#### Usage

    networkSampler$new(type = NA, parameters = NA, nbNodes = NA, directed = FALSE)

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

### `networkSampler$rSamplingMatrix()`

a method for drawing a sampling matrix according to the current sampling
design

#### Usage

    networkSampler$rSamplingMatrix()

------------------------------------------------------------------------

### `networkSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    networkSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
