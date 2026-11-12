# Is the exploration still needed? Benchmark on the blog network

`estimateMissSBM()` fits each number of blocks from a spectral initialization, then polishes
each fit (node swaps) and explores across numbers of blocks (forward/backward passes splitting
and merging blocks, `missSBM_param(iterates = 1)` by default). This benchmark checks whether
these passes are still needed with the variational approximation and the ICL of PR #48, and
whether the PR reduces the local minima.

## Settings

- Instances (`setup.R`): the blog network `frenchblog2007` (196 blogs), fully observed and
  fitted with `node`; 5 samplings with `block-node` as in the JSS analysis (blocks of a
  reference fit with 10 blocks, nodes of the blocks with less than 10% of the blogs observed with
  probability 0.2, the others with 0.8: about 86% of the dyads observed), each fitted with `node`
  and with `block-node`. 11 collections, Q = 1 to 14.
- Settings (`bench.R`): no exploration and no polish (`none`), polish only (`polish`), 1 pass and
  polish (`default`), 2 or 3 passes and polish (`iter2`, `iter3`).
- Versions: `master` (a7cb46d), PR #48 with the ICL of the completed network
  (`PR48-completed`), PR #48 with the ICL built on the observed data (`PR48-observed`, the
  default). `iter2` was only run with `PR48-observed`.
- Metric: for each number of blocks, the gap between the ICL of a setting and the best ICL
  obtained by the same version over all the settings (in practice, `iter3`). It is a gap to the
  best known fit, not to the global optimum. For scale, the ICL is about 6500 on the samplings.

## Reproduce

From the root of the package, with master checked out in a separate worktree:

```sh
git worktree add ../missSBM-master origin/master
export OPENBLAS_NUM_THREADS=1
Rscript inst/simus/exploration-blog/setup.R .
Rscript inst/simus/exploration-blog/bench.R ../missSBM-master master "" 18 none,polish,default,iter3
Rscript inst/simus/exploration-blog/bench.R . PR48-completed completed 18 none,polish,default,iter3
Rscript inst/simus/exploration-blog/bench.R . PR48-observed observed 18
Rscript inst/simus/exploration-blog/summarize.R
```

About 2 minutes per version on 18 cores. The tables are in `results/summary.md`. (`iter2` was
run afterwards, in a separate run, and appended to `results/PR48-observed.rds`.)

## Findings

**The exploration is still needed, with all the versions.** Without it, the ICL is on average
about 400 above the best known, up to 1000 to 1400 for some numbers of blocks, and the curve of
the ICL in Q has 2 to 4 local minima. On the fully observed network, where the three versions
give the same results, 9 blocks are selected instead of 10 (ARI 0.67 with the reference
partition). These local minima come from the initialization and the VEM, not from the missing
data. The polish alone reduces the gap by 10% to 20% only. One pass brings it to about 2 on the
full network, with a single minimum in Q.

**The PR reduces the local minima with `block-node`**, the design whose estimation it changes
(means over the 5 samplings):

| `block-node` | master | PR, completed | PR, observed |
|---|---|---|---|
| ICL gap without exploration | 582 | 476 | 373 |
| Collapsed models without exploration (out of 70) | 60 | 0 | 0 |
| ICL gap after 1 pass (max) | 52 (753) | 18 (190) | 7.6 (95) |
| ICL gap of the selected model after 1 pass | 16.4 | 24.3 | 3.6 |
| Selected Q, no exploration | 9, 11, 12, 14, 11 | 9, 10, 12, 13, 7 | 8, 9, 8, 9, 9 |
| Selected Q, 3 passes | 8, 7, 8, 11, 8 | 8, 11, 11, 8, 9 | 8, 7, 9, 8, 9 |

With the ICL built on the observed data, the selected Q hardly depends on the exploration; on
master it goes from 6 to 14. With `node`, after 1 pass, the gap of the selected model goes from
6.5 (master) to 0.5.

**2 passes give the same selected models as 3**, on the 11 collections, with no gap on the
selected model (gaps remain beyond the selected number of blocks only). With the default
(1 pass), the gap of the selected model is 1.85 on average, and 1 collection out of 11 selects 8
blocks instead of 9.

**Recommendation:** keep `iterates = 1` by default; 2 passes are enough for a reference analysis
such as the JSS one.
