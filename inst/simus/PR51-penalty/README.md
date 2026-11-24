# Simulations for issue #51: the penalty of the ICL counts the observed dyads

The ICL built on the observed data (default since #48) is −2 times the variational bound of
log p(Y^o, R, Z), without the entropy of the blocks, plus a penalty. Its SBM part is the
log-likelihood of the observed dyads, on which the connectivity is estimated. Before #51, the
penalty of the connectivity counted all the possible dyads, log(n(n − 1)/2); it now counts the
observed dyads, as sbm does for a network with missing values. Each connectivity parameter is
penalized −log(rate of observed dyads) less: 1.9 at 15% of observed dyads, 1.0 at 35%.

The risk of a lighter penalty is to select too many blocks, so these simulations target the
regimes where the selection of Q is borderline. Both penalties are computed on the same fits:
the penalty counting all the dyads adds (nbConnectParam + nbCovariates) log(nbDyads /
nbObservedDyads) to the ICL.

## Settings

- Networks: Erdős–Rényi (p = 0.1, true Q = 1); SBM with 2 or 3 balanced blocks, off-diagonal
  connectivity 0.05 and diagonal 0.20/0.25 (Q = 2) or 0.25/0.35 (Q = 3), i.e. close to the
  threshold of detectability once sampled.
- n = 100 and 200; `dyad` sampling at 15% and 35% of observed dyads, `node` sampling at 35% and
  60% of observed nodes.
- 20 replicates per setting (800 networks), Q = 1 to 5 explored with the default controls; the
  exploration is guided by the new ICL.

## Reproduce

From the root of the package:

```sh
OPENBLAS_NUM_THREADS=1 Rscript inst/simus/PR51-penalty/simulate.R 16
Rscript inst/simus/PR51-penalty/summarize.R
```

About 2 minutes on 16 cores. The tables are in `results/summary.md`.

## Findings

- **No structure is found in Erdős–Rényi networks** with either penalty (160 networks, Q = 1
  always selected).
- **The selected Q differs for 7 networks out of 800, all toward the true Q** (SBM with 3 blocks,
  `dyad` sampling at 15% or 35%, or `node` sampling at 35% with n = 100). Elsewhere, both
  penalties select the same Q.
- Near the threshold, the ICL remains conservative with both penalties: at 15% of observed
  dyads, the structure is seldom detected (true Q in 0% to 15% of the cases).
- The campaign of #48 (`inst/simus/PR48-block-designs/`, rerun as `PR51-observed`) is unchanged,
  except for `dyad` sampling in experiment 2: one more network with the true Q (mean selected Q
  2.9 → 3.0, ARI 0.86 → 0.89).

**Choices.** The observed dyads are counted for all the designs, `double-standard` included,
whose connectivity is also informed by the missing dyads through their probability of being
missing (#52), but much less than by the observed ones: the same count for all the designs keeps
the ICLs of the designs comparable on the same data. The covariates of `covar-*` are estimated
on the observed dyads and count them too. The sampling parameters count all the dyads (or all the
nodes), as the sampling matrix is fully observed.
