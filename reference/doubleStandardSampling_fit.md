# Class for fitting a double-standard sampling

A dyad is observed with probability \\p_0\\ if it is 0, \\p_1\\ if it is
1 (the vector of parameters is \\(p_0, p_1)\\). Given the blocks, a dyad
between blocks q and l is missing with probability \\m\_{ql} =
\pi\_{ql}(1 - p_1) + (1 - \pi\_{ql})(1 - p_0)\\: the missing dyads are
integrated out given the blocks, and the posterior probability that a
missing dyad between blocks q and l is an edge is \\\nu\_{ql} =
\pi\_{ql}(1 - p_1) / m\_{ql}\\, a Q x Q matrix. The parameters are
updated jointly with the connectivity of the SBM (see
[SimpleSBM_fit_MNAR](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_MNAR.md)):
this class only stores them, with the terms of the sampling.

With a single block, \\(\pi, p_0, p_1)\\ are not identifiable: the
observed data only inform \\\pi p_1\\ and \\(1 - \pi) p_0\\. The fit is
then constrained to \\p_0 = p_1\\, that is, to dyad sampling, with one
parameter.

## Super classes

[`networkSampling`](https://grosssbm.github.io/missSBM/reference/networkSampling.md)
-\>
[`networkSampling_fit`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.md)
-\>
[`networkSamplingDyads_fit`](https://grosssbm.github.io/missSBM/reference/networkSamplingDyads_fit.md)
-\> `doubleStandardSampling_fit`

## Active bindings

- `df`:

  the number of parameters: 2, or 1 with a single block (p0 = p1)

- `vExpec`:

  variational expectation of the log-likelihood of the sampling, the
  missing dyads being integrated out given the blocks: \\E\_\tau\[\log
  p(R \mid Y^o, Z)\]\\

- `log_lambda`:

  Q x Q matrix, the log-probability that a dyad between blocks q and l
  is missing, the term of the missing dyads in the update of the blocks

- `blockImputation`:

  Q x Q matrix, the posterior probability that a missing dyad between
  blocks q and l is an edge, with the current parameters

## Methods

### Public methods

- [`doubleStandardSampling_fit$new()`](#method-doubleStandardSampling_fit-initialize)

- [`doubleStandardSampling_fit$update_parameters()`](#method-doubleStandardSampling_fit-update_parameters)

- [`doubleStandardSampling_fit$update_imputation()`](#method-doubleStandardSampling_fit-update_imputation)

- [`doubleStandardSampling_fit$reorder()`](#method-doubleStandardSampling_fit-reorder)

- [`doubleStandardSampling_fit$vExpec_at()`](#method-doubleStandardSampling_fit-vExpec_at)

- [`doubleStandardSampling_fit$clone()`](#method-doubleStandardSampling_fit-clone)

Inherited methods

- [`networkSampling$print()`](https://grosssbm.github.io/missSBM/reference/networkSampling.html#method-print)
- [`networkSampling_fit$show()`](https://grosssbm.github.io/missSBM/reference/networkSampling_fit.html#method-show)

------------------------------------------------------------------------

### `doubleStandardSampling_fit$new()`

constructor

#### Usage

    doubleStandardSampling_fit$new(partlyObservedNetwork, blockInit, connectInit)

#### Arguments

- `partlyObservedNetwork`:

  a object with class partlyObservedNetwork representing the observed
  data with possibly missing entries

- `blockInit`:

  n x Q matrix of initial block indicators

- `connectInit`:

  Q x Q matrix of initial connectivity of the SBM

------------------------------------------------------------------------

### `doubleStandardSampling_fit$update_parameters()`

stores the parameters, updated jointly with the connectivity of the SBM
(see
[SimpleSBM_fit_MNAR](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_MNAR.md)),
and the current blocks

#### Usage

    doubleStandardSampling_fit$update_parameters(Z, theta, psi = NULL, ...)

#### Arguments

- `Z`:

  probabilities of block memberships

- `theta`:

  the connectivity of the SBM

- `psi`:

  the updated (p0, p1); `NULL` with a single block, where p0 = p1 is the
  observation rate (dyad sampling)

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `doubleStandardSampling_fit$update_imputation()`

no imputation is needed along the VEM: the missing dyads are integrated
out given the blocks (see the field `blockImputation`)

#### Usage

    doubleStandardSampling_fit$update_imputation(...)

#### Arguments

- `...`:

  use for compatibility

#### Returns

`NULL`

------------------------------------------------------------------------

### `doubleStandardSampling_fit$reorder()`

permute the parameters indexed by the blocks, after the blocks of the
SBM were permuted

#### Usage

    doubleStandardSampling_fit$reorder(o)

#### Arguments

- `o`:

  the permutation applied to the blocks (new block k is former block
  `o[k]`)

------------------------------------------------------------------------

### `doubleStandardSampling_fit$vExpec_at()`

the variational expectation of the log-likelihood of the sampling at
given memberships (used by the ICL, at the MAP memberships)

#### Usage

    doubleStandardSampling_fit$vExpec_at(Z)

#### Arguments

- `Z`:

  an n x Q matrix of memberships

------------------------------------------------------------------------

### `doubleStandardSampling_fit$clone()`

The objects of this class are cloneable with this method.

#### Usage

    doubleStandardSampling_fit$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
