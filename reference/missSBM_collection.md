# An R6 class to represent a collection of SBM fits with missing data

The function
[`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
fits a collection of SBM with missing data for a varying number of
block. These models with class
[`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)
are stored in an instance of an object with class `missSBM_collection`,
described here.

Fields are accessed via active binding and cannot be changed by the
user.

This class comes with a set of R6 methods, some of them being useful for
the user and exported as S3 methods. See the documentation for
[`show()`](https://rdrr.io/r/methods/show.html) and
[`print()`](https://rdrr.io/r/base/print.html)

## Active bindings

- `models`:

  a list of models

- `ICL`:

  the vector of Integrated Classification Criterion (ICL) associated to
  the models in the collection (the smaller, the better)

- `bestModel`:

  the best model according to the ICL, restricted to models without
  collapsed classes when at least one such model is available (see
  `$degenerate`)

- `vBlocks`:

  a vector with the number of blocks

- `occupiedBlocks`:

  a vector with the number of classes actually occupied in each model
  (see
  [missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
  `occupiedBlocks`)

- `degenerate`:

  logical vector, `TRUE` for models with collapsed classes
  (`occupiedBlocks < vBlocks`, see
  [missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
  `repair()`)

- `optimizationSettings`:

  the control list used by estimate()/polish()/explore() when not
  overridden per call (set at construction by
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md))

- `optimizationStatus`:

  a data.frame summarizing the optimization process for all models

## Methods

### Public methods

- [`missSBM_collection$new()`](#method-missSBM_collection-initialize)

- [`missSBM_collection$estimate()`](#method-missSBM_collection-estimate)

- [`missSBM_collection$estimate_chain()`](#method-missSBM_collection-estimate_chain)

- [`missSBM_collection$polish()`](#method-missSBM_collection-polish)

- [`missSBM_collection$explore()`](#method-missSBM_collection-explore)

- [`missSBM_collection$plot()`](#method-missSBM_collection-plot)

- [`missSBM_collection$show()`](#method-missSBM_collection-show)

- [`missSBM_collection$print()`](#method-missSBM_collection-print)

- [`missSBM_collection$clone()`](#method-missSBM_collection-clone)

------------------------------------------------------------------------

### `missSBM_collection$new()`

constructor for networkSampling

#### Usage

    missSBM_collection$new(partlyObservedNet, sampling, clusterInit, control)

#### Arguments

- `partlyObservedNet`:

  An object with class
  [`partlyObservedNetwork`](https://grosssbm.github.io/missSBM/reference/partlyObservedNetwork.md).

- `sampling`:

  The sampling design for the modelling of missing data: MAR designs
  ("dyad", "node") and MNAR designs ("double-standard", "block-dyad",
  "block-node" ,"degree")

- `clusterInit`:

  Initial clustering: a list of vectors, each with size
  `ncol(adjacencyMatrix)`.

- `control`:

  a list of parameters controlling advanced features. Only 'trace' and
  'useCov' are relevant here. See
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
  for details.

------------------------------------------------------------------------

### `missSBM_collection$estimate()`

method to launch the estimation of the collection of models

#### Usage

    missSBM_collection$estimate(control = NULL)

#### Arguments

- `control`:

  optional list of parameters overriding the collection's stored control
  (set at construction by
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md),
  see its details for the full list). Default `NULL` uses the stored
  control as-is.

------------------------------------------------------------------------

### `missSBM_collection$estimate_chain()`

alternative to `estimate()`: fits each model in increasing order of
number of blocks, initializing `vBlocks[k]` by splitting (see
[missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
[`split()`](https://rdrr.io/r/base/split.html)/`candidates_split()`) the
already-converged model at `vBlocks[k-1]` instead of an independent,
cold spectral clustering. Meant to reduce VEM component collapse at
higher numbers of blocks (see `$degenerate`), at the cost of being
sequential in the number of blocks (unlike `estimate()`, which fits
every model in parallel) – can be slower in wall-clock time with many
workers available. Falls back to this slot's own cold-started clustering
(built at construction, same as `estimate()` would use) whenever nothing
is splittable along the chain.

#### Usage

    missSBM_collection$estimate_chain(control = NULL)

#### Arguments

- `control`:

  optional list of parameters overriding the collection's stored control
  (set at construction by
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md),
  see its details for the full list). Default `NULL` uses the stored
  control as-is.

------------------------------------------------------------------------

### `missSBM_collection$polish()`

method to node-swap-polish every model in the collection (see
[missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
`polish()`); fixes individually misclassified nodes at each model's own
number of blocks, unlike `explore()` which searches across blocks.

#### Usage

    missSBM_collection$polish(control = NULL)

#### Arguments

- `control`:

  optional list of parameters overriding the collection's stored control
  (set at construction by
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md),
  see its details for the full list). Default `NULL` uses the stored
  control as-is.

------------------------------------------------------------------------

### `missSBM_collection$explore()`

method for performing exploration of the ICL (split/merge search across
numbers of blocks, see
[missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)'s
`candidates_split()`/`candidates_merge()`). Uses the collection's stored
control by default; `iterates` lets the caller override it for this call
only, without altering the stored control – handy to alternate
`explore()`/`polish()` calls without having to reconstruct a full
control list each time. `iterates <= 0` is a no-op.

#### Usage

    missSBM_collection$explore(control = NULL, iterates = NULL, direction = "both")

#### Arguments

- `control`:

  optional list of parameters overriding the collection's stored control
  (set at construction by
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md),
  see its details for the full list). Default `NULL` uses the stored
  control as-is.

- `iterates`:

  optional integer overriding `control$iterates` for this call only.

- `direction`:

  character ("forward", "backward", "both" or "none") controlling which
  directions are searched. Default "both".

------------------------------------------------------------------------

### `missSBM_collection$plot()`

plot method for missSBM_collection

#### Usage

    missSBM_collection$plot(type = c("icl", "elbo", "monitoring"))

#### Arguments

- `type`:

  the type specifies the field to plot, either "icl", "elbo" or
  "monitoring". Default is "icl"

------------------------------------------------------------------------

### `missSBM_collection$show()`

show method for missSBM_collection

#### Usage

    missSBM_collection$show()

------------------------------------------------------------------------

### `missSBM_collection$print()`

User friendly print method

#### Usage

    missSBM_collection$print()

------------------------------------------------------------------------

### `missSBM_collection$clone()`

The objects of this class are cloneable with this method.

#### Usage

    missSBM_collection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
## Uncomment to set parallel computing with future
## future::plan("multicore", workers = 2)

## Sample 75% of dyads in  French political Blogosphere's network data
adjacencyMatrix <- missSBM::frenchblog2007 %>%
  igraph::delete.vertices(1:100) %>%
  igraph::as_adjacency_matrix() %>%
  missSBM::observeNetwork(sampling = "dyad", parameters = 0.75)
#> Warning: `delete.vertices()` was deprecated in igraph 2.0.0.
#> ℹ Please use `delete_vertices()` instead.
collection <- estimateMissSBM(adjacencyMatrix, 1:5, sampling = "dyad")
#> 
#> 
#>  Adjusting Variational EM for Stochastic Block Model
#> 
#>  Imputation assumes a 'dyad' network-sampling process
#> 
#>  Initialization of 5 model(s). 
#>  Performing VEM inference
#>      Model with 5 blocks.    Model with 2 blocks.    Model with 1 blocks.    Model with 3 blocks.    Model with 4 blocks. Polishing (node-swap)
#> 
#>  Looking for better solutions
#>  Pass 1   Going forward ++++                                                                                                     Pass 1   Going backward ++++                                                                                                    
class(collection)
#> [1] "missSBM_collection" "R6"                
```
