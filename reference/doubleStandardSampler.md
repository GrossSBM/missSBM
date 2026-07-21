# Class for defining a double-standard sampler

Class for defining a double-standard sampler

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampler`](https://grosssbm.github.io/missSBM/reference/networkSampler.md)
-\>
[`dyadSampler`](https://grosssbm.github.io/missSBM/reference/dyadSampler.md)
-\> `doubleStandardSampler`

## Methods

### Public methods

- [`doubleStandardSampler$new()`](#method-doubleStandardSampler-initialize)

- [`doubleStandardSampler$clone()`](#method-doubleStandardSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)
- [`networkSampler$rSamplingMatrix()`](https://grosssbm.github.io/missSBM/reference/networkSampler.html#method-rSamplingMatrix)

------------------------------------------------------------------------

### `doubleStandardSampler$new()`

constructor for networkSampling

#### Usage

    doubleStandardSampler$new(parameters = NA, adjMatrix = NA, directed = FALSE)

#### Arguments

- `parameters`:

  the vector of parameters associated to the sampling at play

- `adjMatrix`:

  matrix of adjacency

- `directed`:

  logical, directed network of not

------------------------------------------------------------------------

### `doubleStandardSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    doubleStandardSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
