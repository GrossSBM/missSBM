# Internal class for a binary SBM fit under the double-standard sampling design (MNAR, depending on the value of each dyad).

The missing dyads are integrated out given the blocks. A dyad between
blocks q and l is missing with probability \\m\_{ql} = \pi\_{ql}(1 -
p_1) + (1 - \pi\_{ql})(1 - p_0)\\, where \\(p_0, p_1)\\ are the
parameters of the sampling: the update of the blocks gets the term
\\\sum\_{j} \sum_l \tau\_{jl} \log m\_{ql}\\ over the missing dyads (i,
j). The missing dyads tie the connectivity to \\(p_0, p_1)\\: the M-step
maximizes over the three jointly, by EM on the missing dyads (the
posterior probability that a missing dyad between blocks q and l is an
edge is \\\nu\_{ql} = \pi\_{ql}(1 - p_1) / m\_{ql}\\), iterated to
convergence, all in Q x Q given the memberships.

It is not designed to be called directly by the user.

## Super classes

[`sbm::SBM`](https://grosssbm.github.io/sbm/reference/SBM.html) -\>
[`sbm::SimpleSBM`](https://grosssbm.github.io/sbm/reference/SimpleSBM.html)
-\>
[`SimpleSBM_fit`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit.md)
-\>
[`SimpleSBM_fit_noCov`](https://grosssbm.github.io/missSBM/reference/SimpleSBM_fit_noCov.md)
-\> `SimpleSBM_fit_MNAR`

## Methods

### Public methods

- [`SimpleSBM_fit_MNAR$update_parameters()`](#method-SimpleSBM_fit_MNAR-update_parameters)

- [`SimpleSBM_fit_MNAR$update_blocks()`](#method-SimpleSBM_fit_MNAR-update_blocks)

- [`SimpleSBM_fit_MNAR$polish_log_tau()`](#method-SimpleSBM_fit_MNAR-polish_log_tau)

- [`SimpleSBM_fit_MNAR$clone()`](#method-SimpleSBM_fit_MNAR-clone)

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

### `SimpleSBM_fit_MNAR$update_parameters()`

update parameters estimation (M-step)

#### Usage

    SimpleSBM_fit_MNAR$update_parameters(psi = NULL, tol = 1e-08, maxit = 10000)

#### Arguments

- `psi`:

  the parameters (p0, p1) of the sampling, starting point of the joint
  M-step. If `NULL`, or with a single block (where the sampling is
  constrained to dyad sampling), the SBM is fitted on the observed
  dyads.

- `tol, maxit`:

  tolerance on the parameters and maximal number of iterations of the EM
  on the missing dyads

#### Returns

invisibly, the updated (p0, p1), or `NULL` if the SBM is fitted on the
observed dyads

------------------------------------------------------------------------

### `SimpleSBM_fit_MNAR$update_blocks()`

update variational estimation of blocks (VE-step)

#### Usage

    SimpleSBM_fit_MNAR$update_blocks(log_lambda = 0)

#### Arguments

- `log_lambda`:

  Q x Q matrix, the log-probability that a dyad between blocks q and l
  is missing (field `log_lambda` of
  [doubleStandardSampling_fit](https://grosssbm.github.io/missSBM/reference/doubleStandardSampling_fit.md))

------------------------------------------------------------------------

### `SimpleSBM_fit_MNAR$polish_log_tau()`

for each node, the log-likelihood of the observed data it would
contribute to each class if hard-assigned there (parameters held fixed),
used to decide node-swap moves in `missSBM_fit$polish()`.

#### Usage

    SimpleSBM_fit_MNAR$polish_log_tau(log_lambda = 0, ...)

#### Arguments

- `log_lambda`:

  Q x Q matrix, the log-probability that a dyad between blocks q and l
  is missing

- `...`:

  additional arguments, not used

#### Returns

an N x Q matrix

------------------------------------------------------------------------

### `SimpleSBM_fit_MNAR$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SimpleSBM_fit_MNAR$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
