# Estimating the degree sampling

In the degree sampling, node i is observed with probability logistic(a + b D_i), D_i its degree.
Its variational bound J_zeta uses a factorized q(Z) q(Y^m), the missing dyads not being
ignorable given the blocks (the sampling of a node depends on its whole row), and the bound of
Jaakkola and Jordan on the sampling term. Three estimators are prototyped in `proto.R`, on top of
the internals of the package (C++ E/M steps, `degreeSampling_fit`), independently of its
estimation; all monitor J_zeta:

- **A**: stage 1, the SBM part on the observed dyads (stops on its bound), then stage 2,
  (nu, a, b, zeta) given the SBM part (stops on J_zeta);
- **B**: the coordinate ascent of J_zeta in all the variables (tau and pi also updated on the
  network completed by nu), from the usual initialization;
- **C**: A, then the coordinate ascent of J_zeta from there (stage 3).

Since 1.2.0, the package implements C (`missSBM_fit`, three stages). The campaign compares the
three prototypes and checks the package against C.

## Settings

- SBM of `model-selection`: Q = 3 blocks of proportions (0.25, 0.35, 0.40), assortative
  connectivity (diagonal 0.6, 0.3, 0.15, off-diagonal 0.05).
- n = 200 and 400; degree sampling with (a, b) = (−0.2, 0), (−1.4, 0.04), (−2.4, 0.08),
  (−4, 0.15) for n = 200 (no, weak, medium, strong effect of the degree; b is divided by 2 for
  n = 400, the degrees doubling), about 41% to 48% of observed nodes in all the settings.
- 20 replicates per setting; the prototypes start from the same spectral clustering, the package
  from its own, at the true Q = 3, without exploration.
- The design is chosen by the ICL among `node`, `block-node` (fitted by the package at Q = 3) and
  `degree` (−2 J_zeta at the MAP memberships + the penalty, as in the package).

## Reproduce

From the root of the package (see `../load-missSBM.R` for `SIMU_LIB`):

```sh
export OPENBLAS_NUM_THREADS=1
SIMU_LIB=<lib-v1.2.0> Rscript inst/simus/degree-estimation/simulate.R 16 200
SIMU_LIB=<lib-v1.2.0> Rscript inst/simus/degree-estimation/simulate.R 16 400
Rscript inst/simus/degree-estimation/summarize.R
```

Less than a minute on 16 cores. The tables are in `results/summary.md`.

## Findings (v1.2.0)

- **C increases the bound in every network**, from the output of A: by 10.0 on average at
  n = 200 (minimum 0.9) and 3.7 at n = 400 (minimum 0.01).
- **C keeps the partition of A**: ARI with A 0.994 to 1.000; the ARI with the truth is the same
  (0.97 to 0.99 at n = 200, 0.998 to 0.999 at n = 400). The feedback loop through the
  imputations of the dyads between unobserved nodes, which makes B fragile, does not appear from
  this starting point.
- **C slightly improves the estimates**: RMSE of the connectivity equal or lower (n = 400,
  strong effect: 0.0053 → 0.0042).
- **With C, the ICL chooses `degree` when it generated the data**: 100% of the cases for a medium
  or strong effect at n = 200 (A: 40% to 50%), 20% for a weak one (A: 0%); at n = 400, 85% to
  100% (A: 75% to 80%) and 60% (A: 40%). Without effect of the degree (b = 0, i.e. node sampling),
  C chooses `degree` in 5% of the cases at n = 200, 0% at n = 400.
- **B is worse**: ARI 0.90 to 0.94 at n = 200, a bound lower than C's by 24 on average (15 at
  n = 400), and `degree` seldom chosen. It reaches a higher bound than C in 1 network out of 80 at
  n = 200, 7 at n = 400: several local maxima.
- **(a, b)** with C, for a medium or strong effect: a within −0.06 and +0.12 of its value, b
  within 3%; with a weak effect, b is underestimated by 10% to 19% for all the estimators
  (b = 0.04 or 0.02: a small parameter, poorly identified).
- **The package matches C**: same ARI on average, same estimates of (a, b) up to 0.02, `degree`
  chosen as often (one network less at n = 200, medium and strong effects). Its ICL is never
  lower than that of C by more than 1, higher in 9 networks out of 160: it starts from its own
  clustering and stops with its own thresholds. It takes 0.2 to 0.35 s per fit.
