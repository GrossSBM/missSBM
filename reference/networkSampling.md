# Definition of R6 Class 'networkSampling'

this virtual class is the mother of all subtypes of networkSampling
(either sampler or fit) It is used to define a sampling model for a
network. It has a rSampling method which takes an adjacency matrix as an
input and send back an object with class partlyObservedNetwork.

## Active bindings

- `type`:

  a character for the type of sampling

- `parameters`:

  the vector of parameters associated with the sampling at play

- `df`:

  the number of entries in the vector of parameters

## Methods

### Public methods

- [`networkSampling$new()`](#method-networkSampling-initialize)

- [`networkSampling$show()`](#method-networkSampling-show)

- [`networkSampling$print()`](#method-networkSampling-print)

- [`networkSampling$clone()`](#method-networkSampling-clone)

------------------------------------------------------------------------

### `networkSampling$new()`

constructor for networkSampling

#### Usage

    networkSampling$new(type = NA, parameters = NA)

#### Arguments

- `type`:

  character for the type of sampling. must be in ("dyad", "covar-dyad",
  "node", "covar-node", "block-node", "block-dyad", "double-standard",
  "degree")

- `parameters`:

  the vector of parameters associated to the sampling at play

------------------------------------------------------------------------

### `networkSampling$show()`

show method

#### Usage

    networkSampling$show(
      type = paste0(private$name, "-model for network sampling\n")
    )

#### Arguments

- `type`:

  character used to specify the type of sampling

------------------------------------------------------------------------

### `networkSampling$print()`

User friendly print method

#### Usage

    networkSampling$print()

------------------------------------------------------------------------

### `networkSampling$clone()`

The objects of this class are cloneable with this method.

#### Usage

    networkSampling$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
