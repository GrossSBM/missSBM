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

- `vExpec`:

  double: variational approximation of the expectation complete
  log-likelihood

- `vExpec_corrected`:

  double: variational approximation of the expectation complete
  log-likelihood with correction to be comparable with MNAR criteria

## Methods

### Public methods

- [`SimpleSBM_fit_withCov$update_parameters()`](#method-SimpleSBM_fit_withCov-update_parameters)

- [`SimpleSBM_fit_withCov$update_blocks()`](#method-SimpleSBM_fit_withCov-update_blocks)

- [`SimpleSBM_fit_withCov$polish_log_tau()`](#method-SimpleSBM_fit_withCov-polish_log_tau)

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

### `SimpleSBM_fit_withCov$polish_log_tau()`

for each node, the complete-data log-likelihood it would contribute to
each class if hard-assigned there (theta/beta/pi held fixed), used to
decide node-swap moves in `missSBM_fit$polish()`.

#### Usage

    SimpleSBM_fit_withCov$polish_log_tau(log_lambda = 0)

#### Arguments

- `log_lambda`:

  additional sampling-design-dependent term, added as-is

#### Returns

an N x Q matrix

------------------------------------------------------------------------

### `SimpleSBM_fit_withCov$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_fit_withCov$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
