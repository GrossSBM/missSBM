# Internal class for a binary SBM fit under MAR sampling designs with covariates.

It is not designed to be called directly by the user.

## Super classes

[`sbm::SBM`](https://grosssbm.github.io/sbm/reference/SBM.html) -\>
[`sbm::SimpleSBM`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html)
-\>
[`SimpleSBM_fit`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.md)
-\> `SimpleSBM_fit_withCov`

## Active bindings

- `imputation`:

  the matrix of imputed values

## Methods

### Public methods

- [`SimpleSBM_fit_withCov$update_parameters()`](#method-SimpleSBM_fit_withCov-update_parameters)

- [`SimpleSBM_fit_withCov$update_blocks()`](#method-SimpleSBM_fit_withCov-update_blocks)

- [`SimpleSBM_fit_withCov$vExpec_at()`](#method-SimpleSBM_fit_withCov-vExpec_at)

- [`SimpleSBM_fit_withCov$clone()`](#method-SimpleSBM_fit_withCov-clone)

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
- [`SimpleSBM_fit$polish_log_tau()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-polish_log_tau)
- [`SimpleSBM_fit$reorder()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-reorder)
- [`SimpleSBM_fit$set_state()`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.html#method-set_state)

------------------------------------------------------------------------

### `SimpleSBM_fit_withCov$update_parameters()`

update parameters estimation (M-step) via Newton-Raphson: the M-step
objective is a weighted logistic regression (concave), so Newton
converges in a handful of iterations – no external optimizer is
required.

#### Usage

    SimpleSBM_fit_withCov$update_parameters(...)

#### Arguments

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `SimpleSBM_fit_withCov$update_blocks()`

update variational estimation of blocks (VE-step)

#### Usage

    SimpleSBM_fit_withCov$update_blocks(...)

#### Arguments

- `...`:

  use for compatibility

------------------------------------------------------------------------

### `SimpleSBM_fit_withCov$vExpec_at()`

variational expectation of the log-likelihood of the observed dyads and
of the blocks at given memberships (used by the ICL, at the MAP
memberships)

#### Usage

    SimpleSBM_fit_withCov$vExpec_at(Z)

#### Arguments

- `Z`:

  an n x Q matrix of memberships

------------------------------------------------------------------------

### `SimpleSBM_fit_withCov$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_fit_withCov$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
