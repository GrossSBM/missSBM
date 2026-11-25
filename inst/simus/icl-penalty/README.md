# The penalty of the ICL counts the observed dyads

The ICL is −2 times the variational bound of log p(Y^o, R, Z) at the MAP memberships, plus a
penalty. Its SBM part is the log-likelihood of the observed dyads, on which the connectivity is
estimated: the penalty of the connectivity counts the observed dyads, as sbm does for a network
with missing values, rather than all the possible dyads, log(n(n − 1)/2). Each connectivity
parameter is then penalized −log(rate of observed dyads) less: 1.9 at 15% of observed dyads, 1.0
at 35%.

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
- 20 replicates per setting (800 networks), Q = 1 to 5 explored with the default controls.

## Reproduce

From the root of the package (see `../load-missSBM.R` for `SIMU_LIB`):

```sh
SIMU_LIB=<lib-v1.2.0> OPENBLAS_NUM_THREADS=1 Rscript inst/simus/icl-penalty/simulate.R 16
Rscript inst/simus/icl-penalty/summarize.R
```

About 1 minute on 16 cores. The tables are in `results/summary.md`.

## Findings (v1.2.0)

- **No structure is found in Erdős–Rényi networks** with either penalty (160 networks, Q = 1
  always selected).
- **The selected Q differs for 13 networks out of 800, all toward the true Q** with the penalty
  on the observed dyads (mostly `dyad` sampling, SBM with 2 or 3 blocks). Elsewhere, both
  penalties select the same Q.
- Near the threshold, the ICL remains conservative with both penalties: at 15% of observed
  dyads, the structure is seldom detected (true Q in 0% to 30% of the cases).

**Choices.** The observed dyads are counted for all the designs, `double-standard` included,
whose connectivity is also informed by the missing dyads through their probability of being
missing, but much less than by the observed ones: the same count for all the designs keeps the
ICLs of the designs comparable on the same data. The covariates of `covar-*` are estimated on the
observed dyads and count them too. The sampling parameters count all the dyads (or all the
nodes), as the sampling matrix is fully observed.
