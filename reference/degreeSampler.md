# Class for defining a degree sampler

Class for defining a degree sampler

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampler`](https://grosssbm.github.io/missSBM/reference/networkSampler.md)
-\>
[`nodeSampler`](https://grosssbm.github.io/missSBM/reference/nodeSampler.md)
-\> `degreeSampler`

## Methods

### Public methods

- [`degreeSampler$new()`](#method-degreeSampler-initialize)

- [`degreeSampler$clone()`](#method-degreeSampler-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-show)
- [`networkSampler$rSamplingMatrix()`](https://grosssbm.github.io/missSBM/reference/networkSampler.html#method-rSamplingMatrix)

------------------------------------------------------------------------

### `degreeSampler$new()`

constructor for networkSampling

#### Usage

    degreeSampler$new(parameters = NA, degrees = NA, directed = FALSE)

#### Arguments

- `parameters`:

  the vector of parameters associated to the sampling at play

- `degrees`:

  vector of nodes' degrees

- `directed`:

  logical, directed network of not

------------------------------------------------------------------------

### `degreeSampler$clone()`

The objects of this class are cloneable with this method.

#### Usage

    degreeSampler$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
