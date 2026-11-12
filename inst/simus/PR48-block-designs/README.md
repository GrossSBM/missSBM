# Simulations for PR #48: block-node and block-dyad

PR #48 fits the SBM part of `block-node` and `block-dyad` on the observed dyads: the missing
dyads, ignorable given the blocks, are integrated out instead of being approximated by a
factorized variational distribution q(Z) q(Y^m). These simulations compare master and the PR
when the number of blocks and the design are unknown. The earlier simulations of the PR fixed
Q to its true value.

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
- **Experiment 2:** data simulated under `node`, `block-node`, `dyad` or `block-dyad` (assortative,
  medium rate; the MAR designs at the mean rate of the block ones). Each network is fitted with
  the MAR design and the block design of its family, and the ICL picks one.

## Reproduce

From the root of the package, with master checked out in a separate worktree:

```sh
git worktree add ../missSBM-master origin/master
export OPENBLAS_NUM_THREADS=1
Rscript inst/simus/PR48-block-designs/simulate.R ../missSBM-master master 16
Rscript inst/simus/PR48-block-designs/simulate.R . PR48 16
Rscript inst/simus/PR48-block-designs/summarize.R master PR48
```

About 3 minutes per version on 16 cores. The results are in `results/`; the tables are in
`results/summary.md`.

## Findings

**At the true number of blocks, the PR is better almost everywhere.** For `block-dyad`
(assortative), the ARI goes from 0.68 to 0.99 at the low rate and from 0.85 to 0.99 at the
medium rate. The RMSE of the connectivity is divided by up to 6 (it is the same at the high
rates). The degenerate fits per collection go from 0.2–1.4 to 0–0.2. Exceptions: the RMSE of
the sampling parameters at the medium rate with `block-dyad` (0.047 → 0.066 assortative,
0.052 → 0.117 disassortative) and at the high rate with `block-node` (0.070 → 0.086
assortative, 0.079 → 0.147 disassortative), and the ARI of `block-node` with disassortative
connectivity at the high rate (0.91 → 0.85).

**With the current ICL, the PR selects too many blocks.** The true Q is selected less often:
for `block-dyad` (assortative, medium rate), 85% on master and 25% with the PR.

**The ICL does not pick the right design, on master as with the PR.** Under `node` sampling,
`block-node` is chosen in 20 cases out of 20. With `dyad`, a MAR design that the PR does not
change, the selected Q is 5.5 on average, for a true Q of 3.

**Cause: the expected log-likelihood of the imputed dyads in the ICL.** The ICL of every design
adds E[log p(Y^m | Z)] (`vExpec_corrected`) on the imputed dyads, without the entropy of their
imputations. As the imputations are computed by the model itself, this term is about
−Σ H(π_{z_i z_j}) over the missing dyads. It rewards purer blocks, i.e. more blocks, with no
data behind them, and the more dyads are missing, the larger the reward.

**An ICL built on the observed data fixes both.** That is, the bound of log p(Y^o, R) without
the entropy of τ: `-2 * (fittedSampling$vExpec + fittedSBM$vExpec) + penalty`. It is computed
here on the same collections, without changing the exploration. With the PR:

- The true Q is selected in 100% of the cases for `block-dyad` (assortative, low and medium
  rates; it was 10% and 25%). For `block-node`, the rate goes from 5% to 85% (medium) and from
  50% to 100% (high).
- The ARI at the selected Q goes from 0.81 to 0.99 (`block-dyad`, low rate) and from 0.75 to
  0.96 (`block-node`, medium rate).
- The true design is chosen in 100% (`node`, it was 0%), 95% (`block-node`), 100% (`dyad`) and
  100% (`block-dyad`) of the cases. The selected Q of `dyad` goes from 5.5 to 3.0.

The weak spots left are `block-node` at the low rate and with disassortative connectivity at the
medium rate, where the true Q is selected in 30% to 50% of the cases: little is observed there.

On master, the ICL built on the observed data has no meaning for the block designs: their SBM
part is fitted on the imputations.

With a single block, this ICL still gives the same value for `node` and `block-node` (resp.
`dyad` and `block-dyad`). For `degree` and `double-standard`, whose sampling part keeps the
imputations, the comparable criterion is the full bound, with the entropy of the imputations.
