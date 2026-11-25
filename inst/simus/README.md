# Simulations of missSBM

Simulation campaigns that check the estimation of missSBM, compare versions and support future
developments. Each campaign has its scripts, a `README.md` (settings, how to reproduce,
findings) and its results. This directory is excluded from the package build (`.Rbuildignore`).

| Directory | Question | Versions |
|---|---|---|
| `model-selection/` | Selection of the number of blocks and of the design by the ICL, with `block-node` and `block-dyad` (missing dyads integrated out given the blocks) and the other designs of each family | v1.1.0, v1.2.0 |
| `exploration/` | Are the polish and the exploration across numbers of blocks needed? Local minima of the ICL in Q, on the French political blogosphere | v1.1.0, v1.2.0 |
| `icl-penalty/` | The penalty of the ICL counts the observed dyads: does it select too many blocks near the threshold of detectability? | v1.2.0 |
| `degree-estimation/` | Estimation of the `degree` sampling: prototypes of three estimators, and the package against the one it implements | v1.2.0 |
| `performance/` | Time and peak memory of the phases of a fit on large sparse networks | v1.1.0, v1.2.0 |

`greed.R` compares missSBM with the package greed on the French political blogosphere.

## Running a campaign

The scripts run from the root of the package and load missSBM through `load-missSBM.R`: the
version installed in the library `SIMU_LIB` if it is set, the sources of the working directory
(`devtools::load_all()`) otherwise. To compare versions, install each one in its own library
with `R CMD INSTALL --preclean` (see `load-missSBM.R`); timings are only meaningful for
installed versions, as `load_all()` compiles the C++ code without optimization. Results are
labelled by the version loaded (`v1.2.0.rds`), and the `summarize.R` scripts compare v1.1.0 and
v1.2.0 by default. Set `OPENBLAS_NUM_THREADS=1`, as the jobs run in parallel.

## Not covered yet, and leads

Gaps of the campaigns, to fill before relying on them for these cases:

- **Covariates**: no campaign on `covar-node`, `covar-dyad`, nor on covariates in the SBM
  (`useCov`). A bug of the M-step with covariates, present since 1.1.0, was only found by the
  replication scripts of the JSS article (`inst/JSS-analyses/`).
- **Directed networks** and **`snowball`**: in no campaign (only in the tests and in the
  numerical checks of the refactorings).

What the campaigns of 1.2.0 suggest:

- **Choice of the design by the ICL**: between designs whose sampling carries the same
  information (`block-node` and `degree` on an assortative network, where the degree depends on
  the block), the ICLs differ by a few units, about the penalty of one parameter (log n), and the
  partitions are the same (`model-selection`). Such a difference is no evidence for one design.
- **`block-node` on the blog network** (`node-vs-block-node-blog.R`, 12 block-node samplings with
  two rates, 0.2 for the small blocks and 0.8 for the others): the ICL prefers `node` in 9 cases
  out of 12 (difference -13 on average, from -69 to +41), with close partitions (ARI with the
  full network 0.68 for `node`, 0.70 for `block-node`) and imputations (AUC 0.83 and 0.84).
  `block-node` has one sampling parameter per block, many more than the two rates of the
  sampling, and the small blocks inform their rate little. A design with fewer rates than blocks
  (rates shared between blocks) would fit such samplings better.
- **Exploration**: one pass (`iterates = 1`, the default) leaves an ICL gap with `node` on the
  block-node samplings of the blog network (8 to 10 blocks), closed by 2 passes (`exploration`);
  on the simulated networks of `model-selection` (Q = 3, n = 200), a second pass changes no
  selected model (ICL changed in 8 fits out of 360, by at most 0.18) and doubles the time
  (`results/summary-v1.2.0_vs_v1.2.0-iter2.md`). Hence the default of 1 pass, and 2 for a
  reference analysis of a real network with many blocks.
- **Local minima**: without exploration, the ICL is 280 to 400 above the best known, in both
  versions: they come from the initialization (spectral clustering) and the VEM, not from the
  missing data. A better initialization is the main lever.
- **Cost**: what remains is the E-step and the M-step on the pattern of the observed dyads,
  nearly full; `degree` costs 3 to 5 times the other designs, in proportion to the number of
  missing dyads (`performance`).
- **`degree`**: with a weak effect of the degree, b is underestimated by 10% to 19% by all the
  estimators (`degree-estimation`): a small parameter, poorly identified.
