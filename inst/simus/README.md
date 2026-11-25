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
