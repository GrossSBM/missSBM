# Time and memory of a fit on large sparse networks

Time and peak memory of the phases of a fit: the network object
(`partlyObservedNetwork$new()`), the initial clustering (`clustering()`), and the fit at the true
number of blocks (`missSBM_fit$doVEM()`, Q = 3), v1.1.0 against v1.2.0. Peaks are the increase of
the R heap during the phase (`gc()` reset before it), i.e. allocations and garbage not yet
collected; memory allocated by C++ code (Armadillo) is not counted.

## Settings

- n = 2000 and 4000, Q = 3 blocks (proportions 0.25, 0.35, 0.40), assortative connectivity
  scaled by 200 / n (mean degree about 26); the input matrix is dense, 30 Mb and 122 Mb.
- Designs: `node` (50% of the nodes), `block-dyad`, `double-standard`, `degree`.
- One network per setting, each in its own process; the VEM with `threshold = 1e-2`,
  `maxIter = 50`, `fixPointIter = 3`.

## Reproduce

From the root of the package, with v1.1.0 and v1.2.0 installed in two libraries, with
`R CMD INSTALL --preclean` (see `../load-missSBM.R`): `devtools::load_all()` compiles the C++
code without optimization, and `R CMD INSTALL` reuses the objects it leaves in `src/`.

```sh
export OPENBLAS_NUM_THREADS=1
SIMU_LIB=<lib-v1.1.0> Rscript inst/simus/performance/measure.R 8
SIMU_LIB=<lib-v1.2.0> Rscript inst/simus/performance/measure.R 8
Rscript inst/simus/performance/summarize.R v1.1.0 v1.2.0
```

The table is in `results/summary.md`.

## Findings (n = 4000)

| | v1.1.0 | v1.2.0 |
|---|---|---|
| `partlyObservedNetwork$new()` | 1.2–2.8 s, +615–716 Mb (5 times the input) | 0.1–0.2 s, +96 Mb (the result) |
| `clustering()` | 0.6–0.7 s, +576–672 Mb | 0.2 s, +91–128 Mb |
| Fit, `node` / `block-dyad` / `double-standard` | 2.9 / 8.7 / 5.7 s, +556 / +821 / +820 Mb | 0.4 / 0.5 / 0.5 s, +10 / +77 / +165 Mb |
| Fit, `degree` | 15.9 s, +1676 Mb | 4.2 s, +452 Mb |
| Total time, `node` / `block-dyad` / `double-standard` / `degree` | 4.7 / 12.1 / 9.0 / 18.0 s | 0.7 / 0.9 / 1.0 / 4.5 s |

- The network object reads the dense input in C++ and builds the sparse patterns of the observed
  dyads, of the missing ones and of the edges directly, without masks, copies or transposes.
- The clustering needs no matrix of distances (farthest pair and distances to the centers in
  C++) nor the indices of the observed dyads.
- The designs that integrate the missing dyads out compute no imputation along the VEM (only on
  demand, at the missing dyads, in C++), and their sampling terms avoid transposing the sparse
  patterns (`crossprod()`).
- The fits of `block-dyad`, `double-standard` and `degree` are not the same algorithm in the two
  versions (missing dyads integrated out, three stages for `degree`): their times compare the
  estimations as a whole. For `degree`, the time and the peak are proportional to the number of
  missing dyads (its imputations), not to n².
- What remains is the E-step and the M-step on the pattern of the observed dyads, nearly full for
  these designs.
