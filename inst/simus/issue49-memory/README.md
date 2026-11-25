# Measurements for issue #49: the remaining dense n x n computations

Time and peak memory of the three phases of a fit: the network object
(`partlyObservedNetwork$new()`), the initial clustering (`clustering()`), and the fit at the true
number of blocks (`missSBM_fit$doVEM()`, Q = 3). Peaks are the increase of the R heap during the
phase (`gc()` reset before it), i.e. allocations and garbage not yet collected; memory allocated
by C++ code (Armadillo) is not counted.

## Settings

- n = 2000 and 4000, Q = 3 blocks (proportions 0.25, 0.35, 0.40), assortative connectivity
  scaled by 200 / n (mean degree about 26); the input matrix is dense, 30 Mb and 122 Mb.
- Designs: `node` (50% of the nodes), `block-dyad`, `double-standard`, `degree`.
- One network per setting, each in its own process.

## Reproduce

From the root of the package:

```sh
OPENBLAS_NUM_THREADS=1 Rscript inst/simus/issue49-memory/baseline.R <label> 8
Rscript inst/simus/issue49-memory/summarize.R   # baseline.rds -> final.rds
```

`results/baseline.rds` is master before #49, `results/final.rds` the end of the branch; the
intermediate files (`point6`, `point5`, `points23`) follow the commits. The table is in
`results/summary.md`.

## Findings (n = 4000)

| | Before | After |
|---|---|---|
| `partlyObservedNetwork$new()` | 1.2–2.8 s, +574–715 Mb (5 times the input) | 0.6–0.7 s, +95 Mb (the result) |
| `clustering()` | 0.6 s, +526–672 Mb | 0.2–0.3 s, +91–128 Mb |
| Fit, `node` / `block-dyad` / `double-standard` | 2.0 / 2.8 / 1.3 s, +553 / +637 / +301 Mb | 1.5 / 1.2 / 1.1 s, +4 / +78 / +110 Mb |
| Total time, `node` / `block-dyad` / `double-standard` | 3.8 / 6.2 / 4.4 s | 2.4 / 2.1 / 2.0 s |
| Fit, `degree` | 10.0 s, +546 Mb | 9.3 s, +523 Mb |

- **Point 6:** the imputations were computed at every E-step for the designs that integrate the
  missing dyads out (13% of the fit time for `node`, 34% for `block-dyad`).
- **Sampling term of `block-dyad` and of `double-standard`:** `t(R) %*% Z` and `t(S) %*% Z`
  formed sparse copies of R and S (about 95 Mb each) at each fixed-point iteration;
  `crossprod()` does not.
- **Point 5:** the masks and indices of the construction (5 times the input) are replaced by a
  reading of the dense matrix in C++, which builds the compressed sparse matrices directly.
- **Points 2 and 3:** no indices of the observed dyads in `imputation()` without covariates; no
  matrix of distances in `kmeans_missSBM()` (farthest pair and distances to the centers in C++).
- **Point 4:** the imputations at the missing dyads in C++, without the n x n matrix.
- **`degree`:** its peak comes from the update of the imputations, which allocates about ten
  vectors of the size of the number of missing dyads (2.5 millions here), and its time from the
  E-step and from `check_boundaries()` on these vectors: proportional to the missing dyads, not
  to n², hence outside the scope of #49.
- **Point 1** (the E-step and the M-step on R, nearly full) is left for #50.
