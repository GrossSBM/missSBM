# Simulations for PR #48: block-node and block-dyad, and the ICL

PR #48 fits the SBM part of `block-node` and `block-dyad` on the observed dyads: the missing
dyads, ignorable given the blocks, are integrated out instead of being approximated by a
factorized variational distribution q(Z) q(Y^m). It also makes the ICL built on the observed
data the default (`missSBM_param(ICL = "observed")`). These simulations compare three versions
when the number of blocks and the design are unknown:

- `master`: before the PR;
- `PR48-completed`: the PR with the ICL of the network completed by the imputations, the
  convention of the previous versions (`ICL = "completed"`);
- `PR48-observed`: the PR with the ICL built on the observed data (the default).

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

From the root of the package, with master checked out in a separate worktree:

```sh
git worktree add ../missSBM-master origin/master
export OPENBLAS_NUM_THREADS=1
Rscript inst/simus/PR48-block-designs/simulate.R ../missSBM-master master 16
SIMU_ICL=completed Rscript inst/simus/PR48-block-designs/simulate.R . PR48-completed 16
SIMU_ICL=observed  Rscript inst/simus/PR48-block-designs/simulate.R . PR48-observed 16
Rscript inst/simus/PR48-block-designs/summarize.R master PR48-completed PR48-observed
```

About 4 minutes per version on 16 cores. The results are in `results/`; the tables are in
`results/summary.md`.

## Findings

**At the true number of blocks, the new variational approximation is better almost everywhere**
(`master` → `PR48-*`, the fit at Q = 3 does not depend on the ICL). For `block-dyad`
(assortative), the ARI goes from 0.68 to 0.99 at the low rate and from 0.85 to 0.99 at the
medium rate; with disassortative connectivity at the low rate, from 0 to 0.54 (0.96 with the
exploration guided by the new ICL). The RMSE of the connectivity is divided by up to 6 (it is
the same at the high rates). The degenerate fits per collection go from 0.2–1.4 to 0–0.3.

**With the ICL of the completed network, the PR selects too many blocks.** For `block-dyad`
(assortative, medium rate), the true Q is selected in 85% of the cases on master and 25% with
the PR. The ICL of the completed network adds −2 E[log p(Y^m | Z)] on the imputed dyads, without
their entropy. As the imputations are computed by the model itself, this term is about
2 Σ H(π_{z_i z_j}) over the missing dyads: it rewards purer blocks, i.e. more blocks, with no
data behind them, and the more dyads are missing, the larger the reward. The new approximation
makes the imputations closer to π_{z_i z_j}, hence the stronger effect.

**The ICL built on the observed data fixes the selection of Q.** The true Q is selected in 95%
to 100% of the cases in 11 settings out of 12 (it was 0% to 100% on master, 5% to 100% with the
completed network), and the ARI at the selected Q is then that at the true Q. The exception is
`block-node` with disassortative connectivity at the low rate (47% of the dyads observed, but
the nodes of the third block observed with probability 0.1): the true Q is selected in 15% of
the cases, 2 blocks are often selected, the ARI is 0.51 (0.32 with the completed network).

**It also fixes the choice of the design**, among the three designs of each family:

- With the completed network (master and PR), `double-standard` is chosen in 20 cases out of 20
  under both `double-standard` and `block-dyad`, and `block-node` in 20 out of 20 under both
  `node` and `degree` (and in 16 or 17 under `block-node`). The right choices for
  `double-standard` and `block-node` only reflect this bias.
- With the observed data, the true design is chosen in 100% (`node`), 95% (`block-node`),
  100% (`dyad`) and 100% (`block-dyad`) of the cases.
- `degree` is never chosen, by any version: `block-node` is chosen instead (18 cases out of 20,
  `node` 2). The assortative network makes the degrees depend on the blocks, so that `block-node`
  is a good approximation of `degree`, with no variational gap (the ICL of `degree` keeps the
  factorized q(Y^m) and the Jaakkola-Jordan bound). The partitions are the same (ARI 0.97).
- `double-standard` is chosen in 60% of the cases (`block-dyad` in the others). Here too, the
  dyads of the dense blocks are the most observed, which `block-dyad` can mimic, and the ICL of
  `double-standard` keeps the gap of the factorized q(Y^m). Issue #52 proposes to integrate the
  missing dyads out given the blocks for this design too, which removes this gap.

With a single block, the ICL built on the observed data still gives the same value for `node`
and `block-node` (resp. `dyad` and `block-dyad`).

## Update for #52: double-standard integrates the missing dyads out

Same campaign, rerun on the branch of #52 (`PR52-observed`, ICL built on the observed data),
which also fixes the order of the sampling parameters of `block-node` and `block-dyad` after the
blocks are reordered. Reproduce with

```sh
SIMU_ICL=observed Rscript inst/simus/PR48-block-designs/simulate.R . PR52-observed 16
Rscript inst/simus/PR48-block-designs/summarize.R master PR48-observed PR52-observed
mv inst/simus/PR48-block-designs/results/summary.md inst/simus/PR48-block-designs/results/summary-PR52.md
git checkout inst/simus/PR48-block-designs/results/summary.md
```

The tables are in `results/summary-PR52.md` (master → PR48-observed → PR52-observed).

- **`double-standard` is now chosen in 100% of the cases** when it generated the data (60% with
  PR #48, `block-dyad` in the others), with Q = 3 selected (3.4 on average before) and an ARI of
  0.99 (0.93). It is never chosen for data simulated under `dyad` or `block-dyad`, whose true
  design is still chosen in 100% of the cases.
- **Experiment 1 is not about `double-standard`**: its changes come from the order of the
  sampling parameters. The RMSE of the sampling parameters at Q = 3 is divided by 2 to 4
  (`block-dyad`: 0.049 → 0.012 at the low rate; `block-node`: 0.095 → 0.050 at the medium rate):
  the parameters returned before were not those of the blocks of the SBM. The price: `polish()`
  no longer adds a permuted sampling term, which acted as a random perturbation and sometimes got
  it out of a local optimum. At the low rate of `block-dyad`, the true Q is selected in 95%
  (assortative) and 85% (disassortative) of the cases, instead of 100% and 95%, and the ARI at
  Q = 3 goes from 0.96 to 0.89 (disassortative). The other settings are unchanged.

## Update for #51: the penalty counts the observed dyads

Rerun as `PR51-observed` (tables in `results/summary-PR51.md`, master → PR52-observed →
PR51-observed). All the cells are unchanged except `dyad` sampling in experiment 2: one more
network with the true Q (mean selected Q 2.9 → 3.0, ARI 0.86 → 0.89). See
`inst/simus/PR51-penalty/` for simulations targeting the regimes where the selection of Q is
borderline.

## Update for #55: degree estimated by the coordinate ascent, ICL at the MAP memberships

Rerun as `PR55` (tables in `results/summary-PR55.md`, PR52-observed → PR51-observed → PR55).
`degree` is now fitted in three stages, ending with the coordinate ascent of its bound, and the
ICL is the bound at the MAP memberships for all the designs.

- **Experiment 2, family of `node`:** `degree` is chosen in 85% of the cases under `degree`
  sampling (0% before); the true design in 95% under `node` (100%) and 60% under `block-node`
  (95%), `degree` being chosen instead: with an assortative network, the observation of a node
  depends on its degree under both designs. Over the three designs, the true one is chosen in
  80% of the cases, against 65% before, when `block-node` was always preferred.
- With the ICL minus the entropy of the soft memberships (the definition before #55),
  `block-node` was recovered in 25% of the cases only: the imputations of `degree`, computed from
  the memberships, make them near-hard, which that definition rewarded.
- **Experiment 1** is unchanged, except `block-node` with disassortative connectivity at the
  low rate: the true Q is selected in 25% of the cases (15%), the ARI at Q = 3 goes from 0.47 to
  0.57 and the RMSE of the sampling parameters from 0.21 to 0.15.
- The family of `dyad` is unchanged.
