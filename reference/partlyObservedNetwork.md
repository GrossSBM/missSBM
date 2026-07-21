# An R6 Class used for internal representation of a partially observed network

This class is not exported to the user

## Active bindings

- `samplingRate`:

  The percentage of observed dyads

- `nbNodes`:

  The number of nodes

- `nbDyads`:

  The number of dyads

- `is_directed`:

  logical indicating if the network is directed or not

- `networkData`:

  The adjacency matrix of the network

- `covarArray`:

  the array of covariates

- `covarMatrix`:

  the matrix of covariates

- `samplingMatrix`:

  matrix of observed and non-observed edges

- `samplingMatrixBar`:

  matrix of observed and non-observed edges

- `observedNodes`:

  a vector of observed and non-observed nodes (observed means at least
  one non NA value)

## Methods

### Public methods

- [`partlyObservedNetwork$new()`](#method-partlyObservedNetwork-initialize)

- [`partlyObservedNetwork$clustering()`](#method-partlyObservedNetwork-clustering)

- [`partlyObservedNetwork$imputation()`](#method-partlyObservedNetwork-imputation)

- [`partlyObservedNetwork$clone()`](#method-partlyObservedNetwork-clone)

------------------------------------------------------------------------

### `partlyObservedNetwork$new()`

constructor

#### Usage

    partlyObservedNetwork$new(
      adjacencyMatrix,
      covariates = list(),
      similarity = l1_similarity
    )

#### Arguments

- `adjacencyMatrix`:

  The adjacency matrix of the network

- `covariates`:

  A list with M entries (the M covariates), each of whom being either a
  size-N vector or N x N matrix.

- `similarity`:

  An R x R -\> R function to compute similarities between node
  covariates. Default is `l1_similarity`, that is, -abs(x-y).

------------------------------------------------------------------------

### `partlyObservedNetwork$clustering()`

method to cluster network data with missing value

#### Usage

    partlyObservedNetwork$clustering(
      vBlocks,
      imputation = ifelse(is.null(private$phi), "median", "average")
    )

#### Arguments

- `vBlocks`:

  The vector of number of blocks considered in the collection.

- `imputation`:

  character indicating the type of imputation among "median", "average"

------------------------------------------------------------------------

### `partlyObservedNetwork$imputation()`

basic imputation from existing clustering

#### Usage

    partlyObservedNetwork$imputation(type = c("median", "average", "zero"))

#### Arguments

- `type`:

  a character, the type of imputation. Either "median" or "average"

------------------------------------------------------------------------

### `partlyObservedNetwork$clone()`

The objects of this class are cloneable with this method.

#### Usage

    partlyObservedNetwork$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
