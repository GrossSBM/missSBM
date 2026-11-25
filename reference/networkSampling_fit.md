# Virtual class of the fits of a sampling design

The mother of
[`networkSamplingDyads_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.md)
and
[`networkSamplingNodes_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingNodes_fit.md),
with the defaults of a MAR design: no parameter indexed by the blocks,
nothing to update along the VEM.

## Super class

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\> `networkSampling_fit`

## Active bindings

- `penalty`:

  double, value of the penalty term in ICL: the parameters of the
  sampling count all the sampling units (dyads or nodes), the sampling
  matrix being fully observed

- `log_lambda`:

  the term of the sampling in the update of the blocks: 0 when it does
  not depend on them

## Methods

### Public methods

- [`networkSampling_fit$show()`](#method-networkSampling_fit-show)

- [`networkSampling_fit$update_parameters()`](#method-networkSampling_fit-update_parameters)

- [`networkSampling_fit$update_imputation()`](#method-networkSampling_fit-update_imputation)

- [`networkSampling_fit$reorder()`](#method-networkSampling_fit-reorder)

- [`networkSampling_fit$vExpec_at()`](#method-networkSampling_fit-vExpec_at)

- [`networkSampling_fit$clone()`](#method-networkSampling_fit-clone)

Inherited methods

- [`networkSampling$initialize()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-initialize)
- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)

------------------------------------------------------------------------

### `networkSampling_fit$show()`

show method

#### Usage

    networkSampling_fit$show()

------------------------------------------------------------------------

### `networkSampling_fit$update_parameters()`

update the parameters along the VEM: nothing to do by default, the
parameters of a MAR design not depending on the SBM

#### Usage

    networkSampling_fit$update_parameters(...)

#### Arguments

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `networkSampling_fit$update_imputation()`

a method to update the imputation of the missing entries.

#### Usage

    networkSampling_fit$update_imputation(nu)

#### Arguments

- `nu`:

  the matrix of (uncorrected) imputation for missing entries

------------------------------------------------------------------------

### `networkSampling_fit$reorder()`

permute the parameters indexed by the blocks, after the blocks of the
SBM were permuted. By default, nothing to do (no parameter indexed by
the blocks)

#### Usage

    networkSampling_fit$reorder(o)

#### Arguments

- `o`:

  the permutation applied to the blocks (new block k is former block
  `o[k]`)

------------------------------------------------------------------------

### `networkSampling_fit$vExpec_at()`

the variational expectation of the log-likelihood of the sampling at
given memberships (used by the ICL, at the MAP memberships)

#### Usage

    networkSampling_fit$vExpec_at(Z)

#### Arguments

- `Z`:

  an n x Q matrix of memberships

------------------------------------------------------------------------

### `networkSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    networkSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
