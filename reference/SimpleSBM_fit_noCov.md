# Internal class for a binary SBM fit on the observed dyads, without covariates: MAR sampling designs, 'block-node' and 'block-dyad' (with their sampling term in the update of the blocks) and 'degree'.

It is not designed to be called directly by the user.

## Super classes

[`sbm::SBM`](https://grosssbm.github.io/sbm/reference/SBM.html) -\>
[`sbm::SimpleSBM`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html)
-\>
[`SimpleSBM_fit`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.md)
-\> `SimpleSBM_fit_noCov`

## Active bindings

- `imputation`:

  the matrix of imputed values

- `vExpec`:

  double: variational approximation of the expectation complete
  log-likelihood

## Methods

### Public methods

- [`SimpleSBM_fit_noCov$update_parameters()`](#method-SimpleSBM_fit_noCov-update_parameters)

- [`SimpleSBM_fit_noCov$update_blocks()`](#method-SimpleSBM_fit_noCov-update_blocks)

- [`SimpleSBM_fit_noCov$polish_log_tau()`](#method-SimpleSBM_fit_noCov-polish_log_tau)

- [`SimpleSBM_fit_noCov$vExpec_imputed()`](#method-SimpleSBM_fit_noCov-vExpec_imputed)

- [`SimpleSBM_fit_noCov$clone()`](#method-SimpleSBM_fit_noCov-clone)

Inherited methods

- [`sbm::SBM$print()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-print)
- [`sbm::SBM$rNetwork()`](https://grosssbm.github.io/sbm/reference/SBM.html#method-rNetwork)
- [`sbm::SimpleSBM$plot()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-plot)
- [`sbm::SimpleSBM$predict()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-predict)
- [`sbm::SimpleSBM$rEdges()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rEdges)
- [`sbm::SimpleSBM$rMemberships()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-rMemberships)
- [`sbm::SimpleSBM$show()`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html#method-show)
- [`SimpleSBM_fit$doVEM()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-doVEM)
- [`SimpleSBM_fit$get_state()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-get_state)
- [`SimpleSBM_fit$initialize()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-initialize)
- [`SimpleSBM_fit$reorder()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-reorder)
- [`SimpleSBM_fit$set_state()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-set_state)

------------------------------------------------------------------------

### `SimpleSBM_fit_noCov$update_parameters()`

update parameters estimation (M-step)

#### Usage

    SimpleSBM_fit_noCov$update_parameters(...)

#### Arguments

- `...`:

  additional arguments, only required for MNAR cases

------------------------------------------------------------------------

### `SimpleSBM_fit_noCov$update_blocks()`

update variational estimation of blocks (VE-step)

#### Usage

    SimpleSBM_fit_noCov$update_blocks(log_lambda = 0, ...)

#### Arguments

- `log_lambda`:

  term of the sampling design depending on the blocks ('block-node' and
  'block-dyad'), added to the update of the blocks; 0 for the other
  designs

- `...`:

  additional arguments, only required for MNAR cases

------------------------------------------------------------------------

### `SimpleSBM_fit_noCov$polish_log_tau()`

for each node, the complete-data log-likelihood it would contribute to
each class if hard-assigned there (theta/pi held fixed), used to decide
node-swap moves in `missSBM_fit$polish()`.

#### Usage

    SimpleSBM_fit_noCov$polish_log_tau(log_lambda = 0)

#### Arguments

- `log_lambda`:

  additional sampling-design-dependent term, added as-is

#### Returns

an N x Q matrix

------------------------------------------------------------------------

### `SimpleSBM_fit_noCov$vExpec_imputed()`

variational expectation of the log-likelihood of the missing dyads,
given their imputations (used by the ICL of 'degree', see
[missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md))

#### Usage

    SimpleSBM_fit_noCov$vExpec_imputed(nu)

#### Arguments

- `nu`:

  the imputations of the missing dyads (sparse matrix with the pattern
  of the missing dyads)

------------------------------------------------------------------------

### `SimpleSBM_fit_noCov$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_fit_noCov$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
