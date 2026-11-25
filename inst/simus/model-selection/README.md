# Model selection: number of blocks and design

Since 1.2.0, the designs that depend on the blocks (`block-node`, `block-dyad`) and
`double-standard` integrate the missing dyads out given the blocks, instead of a factorized
variational distribution q(Z) q(Y^m); `degree` maximizes the variational bound of its model; the
ICL is built on the observed data, evaluated at the MAP memberships, with a penalty counting the
observed dyads. These simulations check the selection of the number of blocks and of the design
by the ICL, v1.1.0 against v1.2.0.

## Settings

- SBM with n = 200 nodes and Q = 3 blocks of proportions (0.25, 0.35, 0.40). Connectivity is
  either assortative (diagonal 0.6, 0.3, 0.15, off-diagonal 0.05) or disassortative (see
  `simulate.R`).
- Three observation rates for each design (low, medium, high), given by the block parameters
  ρ (`block-node`) and ψ (`block-dyad`) in `simulate.R`.
- 20 replicates per setting, Q = 1 to 6 explored with the default controls. The data are the
  same for both versions.
- **Experiment 1:** data simulated and fitted with `block-node` or `block-dyad`, 2 connectivities
  x 3 rates.
- **Experiment 2:** data simulated under one of the designs of a family, `node`, `block-node`,
  `degree` or `dyad`, `block-dyad`, `double-standard` (assortative connectivity, medium rate;
  the MAR designs at the mean rate of the block ones; `degree` with ψ = (−2.4, 0.08) and
  `double-standard` with (p0, p1) = (0.3, 0.7), close to the rates of the other designs of the
  family). Each network is fitted with the three designs of its family, and the ICL picks one.

## Reproduce

From the root of the package, with v1.1.0 and v1.2.0 installed in two libraries (see
`../load-missSBM.R`):

```sh
export OPENBLAS_NUM_THREADS=1
SIMU_LIB=<lib-v1.1.0> Rscript inst/simus/model-selection/simulate.R 16
SIMU_LIB=<lib-v1.2.0> Rscript inst/simus/model-selection/simulate.R 16
Rscript inst/simus/model-selection/summarize.R v1.1.0 v1.2.0
```

About 4 minutes (v1.1.0) and 2.5 minutes (v1.2.0) on 16 cores. The tables are in
`results/summary.md`.

## Findings

**At the true number of blocks, the estimation of the block designs is much better.** For
`block-dyad` at the low rate, the ARI at Q = 3 goes from 0.68 to 0.97 (assortative) and from 0 to
0.89 (disassortative); for `block-node` at the low and medium rates, from 0.19–0.78 to
0.57–0.99. At these rates, the RMSE of the connectivity and that of the sampling parameters are
divided by 2.5 to 5.5; at the high rates, both versions are close. Collapsed fits almost disappear (0 to 0.3 per
collection, against 0.2 to 1.4).

**The true number of blocks is selected**: in 95% to 100% of the cases in 10 settings out of 12
(85% for `block-dyad`, disassortative, low rate), against 0% to 100% with v1.1.0, which selected
too many blocks for `block-node` (3.6 to 4.5 on average at the low and medium rates). The exception is
`block-node` with disassortative connectivity at the low rate, where the nodes of the third block
are observed with probability 0.1: Q = 3 in 25% of the cases (2 blocks are often selected), ARI
0.53.

**The ICL picks the design that generated the data**, among the three designs of each family:
`node` 95%, `block-node` 60%, `degree` 90%, `dyad` 100%, `block-dyad` 100%, `double-standard`
100%, against 0% for `node`, `degree` and `block-dyad` with v1.1.0, whose ICL, computed on the
network completed by the imputations, favoured the designs with more imputations
(`block-node`, `double-standard`). Under `block-node`, `degree` is chosen in the other 40%: with an
assortative network, the observation of a node depends on its degree under both designs, and the
three fits give the same partition (ARI 0.97). The selected Q is 3 with the true design in all the
cases (3.0 to 5.2 with v1.1.0).

**A second pass of the exploration changes nothing here** (`SIMU_ITERATES=2`, tables in
`results/summary-v1.2.0_vs_v1.2.0-iter2.md`): the same selected models in all the cases, the ICL
changed in 8 fits out of 360 by at most 0.18, for twice the time.

```sh
SIMU_LIB=<lib-v1.2.0> SIMU_ITERATES=2 Rscript inst/simus/model-selection/simulate.R 16
Rscript inst/simus/model-selection/summarize.R v1.2.0 v1.2.0-iter2
```
