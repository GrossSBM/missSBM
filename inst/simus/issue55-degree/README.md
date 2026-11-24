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
`degreeSampling_fit`), without changing it. All monitor **J_zeta**. (The package currently uses
for the sampling term log logistic(a + b E[D_i]), the expected degree plugged in: as
log logistic is concave, this is above E[log logistic(a + b D_i)], hence not a lower bound.)

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

- **C increases the bound in every network**, from the output of A: by 9.9 on average at
  n = 200 (minimum 0.7) and 3.7 at n = 400 (minimum 0.01). Its decreases along the iterations
  are numerical (relative size below 10^-6).
- **C keeps the partition of A**: ARI with A 0.994 to 1.000; the ARI with the truth is the same
  (0.97 to 0.99 at n = 200, 0.998 to 0.999 at n = 400). The feedback loop through the
  imputations of the dyads between unobserved nodes, which made the thesis VEM fragile, does not
  appear from this starting point.
- **C slightly improves the estimates**: RMSE of the connectivity equal or lower (n = 400,
  strong effect: 0.0053 → 0.0042), a and b slightly less biased.
- **With C, the ICL chooses `degree` when it generated the data**: 75% to 100% of the cases at
  n = 200 (A: 0% to 20%), 60% to 100% at n = 400 (A: 20% to 70%), depending on the effect.
  Without effect of the degree (b = 0, i.e. node sampling), C chooses `degree` in 20% of the
  cases at n = 200 (A: 0%), 0% at n = 400.
- **B, the thesis from the usual initialization, is worse**: ARI 0.90 to 0.94 at n = 200, a
  bound lower than C's by 25 on average (15 at n = 400), and `degree` seldom chosen. It reaches a
  higher bound than C in 1 network out of 80 at n = 200, 6 at n = 400: several local maxima.
- **Cost**: C takes 3 to 4 times A (0.2 s at n = 200, 0.5 s at n = 400, Q = 3).
- **All three estimators bias (a, b) when the effect of the degree is strong**: a overestimated
  by about 0.5 and b underestimated by about 15%, for A, B and C alike. This comes from the model
  or its variational approximation (factorized q(Y^m), Jaakkola–Jordan bound), not from the
  algorithm: to be investigated separately.

## Consequences for #55

The simulations support option C. For the implementation:

- the bound of `degree` must be J_zeta: replace the plug-in sampling term of
  `degreeSampling_fit$vExpec` by the Jaakkola–Jordan bound, in the monitored objective and in the
  ICL (the ICL of `degree` changes even before stage 3);
- stage 3 needs the SBM part updated on the network completed by nu, on the sparse pattern of
  the missing dyads;
- the comparisons here are at the true Q, without exploration: the cost of the three stages for
  the candidates of the exploration, and the choice of Q, remain to be checked in the package.
