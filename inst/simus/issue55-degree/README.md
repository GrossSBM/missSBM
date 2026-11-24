# Simulations for issue #55: estimating the degree sampling

The model of the degree sampling (node i observed with probability logistic(a + b D_i), D_i its
degree) and its variational bound J_zeta (factorized q(Z) q(Y^m), Jaakkola–Jordan bound on the
sampling; Formularius, `eq-bound-degree`) are those of the thesis. The algorithm of the package is
not (Formularius, `sec-degree-estimation`): since missSBM 1.0.0, the SBM part is the VEM of a MAR
design on the observed dyads, which ignores the sampling, and the sampling part is fitted given
it. These simulations compare three estimators before changing the package:

- **A**: stage 1, the SBM part on the observed dyads (stops on its bound), then stage 2,
  (nu, a, b, zeta) given the SBM part (stops on J_zeta). This is the current estimator, with
  justified stopping rules.
- **B**: the coordinate ascent of J_zeta in all the variables (tau and pi also updated on the
  network completed by nu), from the usual initialization: the thesis.
- **C**: A, then the coordinate ascent of J_zeta from there (stage 3). The target of #55.

The three are prototyped in `proto.R` on top of the internals of the package (C++ E/M steps,
`degreeSampling_fit`). All monitor **J_zeta**.

These simulations found two bugs in the bound of the sampling, fixed in `e6a3d08` (the results
below are computed after the fix; those before are in the history, commit `edaa6b9`):

- zeta_i, the parameter of the Jaakkola–Jordan bound, at its optimum sqrt(E[(a + b D_i)^2]),
  often above 1, was truncated to [eps, 1 − eps] like a probability (since 2018). This biased
  (a, b) for all the estimators: with a strong effect of the degree, a was overestimated by 0.6
  and b underestimated by 16% at n = 200; after the fix, 0.03 and 0.6%.
- `degreeSampling_fit$vExpec` plugged the expected degree into log logistic(a + b D_i), which is
  above E[log logistic(a + b D_i)] by Jensen's inequality: not a lower bound. It is now the
  Jaakkola–Jordan bound, used in the ICL.

## Settings

- SBM of the campaign of #48: Q = 3 blocks of proportions (0.25, 0.35, 0.40), assortative
  connectivity (diagonal 0.6, 0.3, 0.15, off-diagonal 0.05).
- n = 200 and 400; degree sampling with (a, b) = (−0.2, 0), (−1.4, 0.04), (−2.4, 0.08),
  (−4, 0.15) for n = 200 (no, weak, medium, strong effect of the degree; b is divided by 2 for
  n = 400, the degrees doubling), about 41% to 48% of observed nodes in all the settings.
- 20 replicates per setting; the three estimators start from the same spectral clustering, at
  the true Q = 3, without exploration.
- The design is chosen by the ICL among `node`, `block-node` (fitted by the package at Q = 3) and
  `degree` (−2 (J_zeta − H(tau)) + the penalty of the package).

## Reproduce

From the root of the package:

```sh
OPENBLAS_NUM_THREADS=1 Rscript inst/simus/issue55-degree/simulate.R 16 200
OPENBLAS_NUM_THREADS=1 Rscript inst/simus/issue55-degree/simulate.R 16 400
Rscript inst/simus/issue55-degree/summarize.R
```

Less than a minute on 16 cores. The tables are in `results/summary.md`.

## Findings

- **C increases the bound in every network**, from the output of A: by 10.0 on average at
  n = 200 (minimum 0.9) and 3.7 at n = 400 (minimum 0.01). Its decreases along the iterations
  are numerical (relative size below 10^-6).
- **C keeps the partition of A**: ARI with A 0.994 to 1.000; the ARI with the truth is the same
  (0.97 to 0.99 at n = 200, 0.998 to 0.999 at n = 400). The feedback loop through the
  imputations of the dyads between unobserved nodes, which made the thesis VEM fragile, does not
  appear from this starting point.
- **C slightly improves the estimates**: RMSE of the connectivity equal or lower (n = 400,
  strong effect: 0.0053 → 0.0042), a and b slightly less biased at n = 400.
- **With C, the ICL chooses `degree` when it generated the data** (ICL at the MAP memberships,
  as in the package since #55): 100% of the cases for a medium or strong effect at n = 200 (A:
  40% to 50%), 20% for a weak one (A: 0%); at n = 400, 85% to 100% (A: 75% to 80%) and 60% (A:
  40%). Without effect of the degree (b = 0, i.e. node sampling), C chooses `degree` in 5% of the
  cases at n = 200 (A: 0%), 0% at n = 400. (With the ICL minus the entropy of the soft
  memberships, used first, C chose `degree` in 75% of the cases with a weak effect and 20%
  without effect at n = 200: the imputations make the memberships near-hard, which that
  definition rewarded.)
- **B, the thesis from the usual initialization, is worse**: ARI 0.90 to 0.94 at n = 200, a
  bound lower than C's by 25 on average (15 at n = 400), and `degree` seldom chosen. It reaches a
  higher bound than C in 1 network out of 80 at n = 200, 6 at n = 400: several local maxima.
- **Cost**: C takes 3 to 4 times A (0.2 s at n = 200, 0.5 s at n = 400, Q = 3).
- **(a, b)** are now estimated without bias for a strong effect at n = 200 (A: a +0.03,
  b −0.6%); at n = 400, a remains overestimated by 0.1 to 0.2 and b underestimated by 2% to 5%;
  with a weak effect, b is underestimated by 10% to 19% for all the estimators (b = 0.04 or 0.02:
  a small parameter, poorly identified).

## Consequences for #55

The simulations support option C. For the implementation:

- stage 3 needs the SBM part updated on the network completed by nu, on the sparse pattern of
  the missing dyads;
- the comparisons here are at the true Q, without exploration: the cost of the three stages for
  the candidates of the exploration, and the choice of Q, remain to be checked in the package.
