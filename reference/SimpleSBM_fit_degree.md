# Internal class for a binary SBM fit under the degree sampling design (MNAR, depending on the degrees).

Without imputations (`nu = NULL`), the SBM is fitted on the observed
dyads, as
[`SimpleSBM_fit_noCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.md):
the first stage of the estimation of this design. With the imputations
of the sampling, the blocks and the connectivity are updated on the
network completed by them: the coordinate ascent of the variational
bound of the degree sampling (factorized q(Z) q(Y^m)), its third stage
(see
[missSBM_fit](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)).
Only the pattern of the missing dyads is used for the imputations.

It is not designed to be called directly by the user.

## Super classes

[`sbm::SBM`](https://grosssbm.github.io/sbm/reference/SBM.html) -\>
[`sbm::SimpleSBM`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html)
-\>
[`SimpleSBM_fit`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.md)
-\>
[`SimpleSBM_fit_noCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.md)
-\> `SimpleSBM_fit_degree`

## Methods

### Public methods

- [`SimpleSBM_fit_degree$update_parameters()`](#method-SimpleSBM_fit_degree-update_parameters)

- [`SimpleSBM_fit_degree$update_blocks()`](#method-SimpleSBM_fit_degree-update_blocks)

- [`SimpleSBM_fit_degree$polish_log_tau()`](#method-SimpleSBM_fit_degree-polish_log_tau)

- [`SimpleSBM_fit_degree$clone()`](#method-SimpleSBM_fit_degree-clone)

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
- [`SimpleSBM_fit_noCov$vExpec_at()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.html#method-vExpec_at)
- [`SimpleSBM_fit_noCov$vExpec_imputed()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.html#method-vExpec_imputed)

------------------------------------------------------------------------

### `SimpleSBM_fit_degree$update_parameters()`

update parameters estimation (M-step)

#### Usage

    SimpleSBM_fit_degree$update_parameters(nu = NULL, ...)

#### Arguments

- `nu`:

  the imputations of the missing dyads (sparse matrix with the pattern
  of the missing dyads), or `NULL` to fit the SBM on the observed dyads

- `...`:

  additional arguments, not used

------------------------------------------------------------------------

### `SimpleSBM_fit_degree$update_blocks()`

update variational estimation of blocks (VE-step)

#### Usage

    SimpleSBM_fit_degree$update_blocks(log_lambda = 0, nu = NULL, ...)

#### Arguments

- `log_lambda`:

  not used (0 for the degree sampling)

- `nu`:

  the imputations of the missing dyads, or `NULL` to use the observed
  dyads only

- `...`:

  additional arguments, not used

------------------------------------------------------------------------

### `SimpleSBM_fit_degree$polish_log_tau()`

for each node, the expected complete log-likelihood it would contribute
to each class if hard-assigned there (parameters and imputations held
fixed), used to decide node-swap moves in `missSBM_fit$polish()`.

#### Usage

    SimpleSBM_fit_degree$polish_log_tau(log_lambda = 0, nu = NULL, ...)

#### Arguments

- `log_lambda`:

  not used (0 for the degree sampling)

- `nu`:

  the imputations of the missing dyads, or `NULL` to use the observed
  dyads only

- `...`:

  additional arguments, not used

#### Returns

an N x Q matrix

------------------------------------------------------------------------

### `SimpleSBM_fit_degree$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_fit_degree$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
