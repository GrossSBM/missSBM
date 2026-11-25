# Is the exploration needed? Benchmark on the blog network

`estimateMissSBM()` fits each number of blocks from a spectral initialization, then polishes
each fit (node swaps) and explores across numbers of blocks (forward/backward passes splitting
and merging blocks, `missSBM_param(iterates = 1)` by default). This benchmark checks whether
these passes are needed, how many, and how often the ICL has local minima in Q, v1.1.0 against
v1.2.0.

## Settings

- Instances (`setup.R`, drawn once and kept fixed, `results/instances.rds`): the blog network
  `frenchblog2007` (196 blogs), fully observed and fitted with `node`; 5 samplings with
  `block-node` as in the JSS analysis (blocks of a reference fit with 10 blocks, nodes of the
  blocks with less than 10% of the blogs observed with probability 0.2, the others with 0.8:
  about 86% of the dyads observed), each fitted with `node` and with `block-node`.
  11 collections, Q = 1 to 14.
- Settings (`bench.R`): no exploration and no polish (`none`), polish only (`polish`), 1 pass and
  polish (`default`), 2 or 3 passes and polish (`iter2`, `iter3`).
- Metric: for each number of blocks, the gap between the ICL of a setting and the best ICL
  obtained by the same version over all the settings (in practice, `iter3`). It is a gap to the
  best known fit, not to the global optimum. For scale, the ICL is about 6500 on the samplings.

## Reproduce

From the root of the package, with v1.1.0 and v1.2.0 installed in two libraries (see
`../load-missSBM.R`):

```sh
export OPENBLAS_NUM_THREADS=1
SIMU_LIB=<lib-v1.1.0> Rscript inst/simus/exploration/bench.R 16
SIMU_LIB=<lib-v1.2.0> Rscript inst/simus/exploration/bench.R 16
Rscript inst/simus/exploration/summarize.R v1.1.0 v1.2.0
```

About 2 minutes per version on 16 cores. The tables are in `results/summary.md`.

## Findings

**The exploration is needed, with both versions.** Without it, the ICL is on average 280 to 580
above the best known, up to 900 to 1400 for some numbers of blocks, and its curve in Q has 2 to 4
local minima. On the fully observed network, where both versions give the same results, 9 blocks
are selected instead of 10 (ARI 0.67 with the reference partition). These local minima come from
the initialization and the VEM, not from the missing data. The polish alone reduces the gap by a
few percent only; one pass brings it to about 2 on the full network, with a single minimum in Q.

**With `block-node`, the design whose estimation changed most, v1.2.0 has fewer local minima**
(means over the 5 samplings):

| `block-node` | v1.1.0 | v1.2.0 |
|---|---|---|
| ICL gap without exploration | 582 | 366 |
| Collapsed models without exploration (out of 70) | 60 | 0 |
| ICL gap after 1 pass (max) | 52 (753) | 8.4 (76) |
| ICL gap of the selected model after 1 pass | 16.4 | 1.0 |
| Selected Q, no exploration | 9, 11, 12, 14, 11 | 8, 9, 8, 8, 9 |
| Selected Q, 3 passes | 8, 7, 8, 11, 8 | 8, 8, 8, 6, 9 |

With v1.2.0, the selected Q hardly depends on the exploration; with v1.1.0 it goes from 6 to 14.

**With `node` on the samplings, one pass leaves a larger gap with v1.2.0** (12 on the selected
model, 6.5 with v1.1.0), and the selected Q changes with a second pass in 1 sampling out of 5;
2 passes bring the gap of the selected model to about 3, 3 passes to 0. v1.2.0 selects 8 or 9
blocks, against 9 or 10 with v1.1.0 (ICL built on the observed data rather than on the network
completed by the imputations).

**Cost:** v1.2.0 is about twice as fast (one pass: 8 to 11 s per collection, against 18 to 26 s;
16 jobs in parallel).

**Recommendation:** keep `iterates = 1` by default; 2 passes for a reference analysis such as the
JSS one.
