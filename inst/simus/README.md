# Simulations of missSBM

Simulation campaigns that support the changes of the estimation, each with its scripts, a
`README.md` (settings, how to reproduce, findings) and its results. This directory is excluded
from the package build (`.Rbuildignore`).

| Directory | Supports | Question |
|---|---|---|
| `PR48-block-designs/` | #48, rerun for #52, #51, #55 | `block-node` and `block-dyad` with the missing dyads integrated out; the ICL built on the observed data; the choice of the design and of the number of blocks |
| `exploration-blog/` | #48 | the exploration of the number of blocks on the French political blogosphere |
| `PR51-penalty/` | #51 | the penalty of the ICL counting the observed dyads, near the threshold of detectability |
| `issue55-degree/` | #55 | the estimation of the `degree` sampling: two stages, the thesis, two stages then the coordinate ascent |
| `issue49-memory/` | #49 | time and peak memory of the network object, the initial clustering and the fit |

`greed.R` compares missSBM with the package greed on the French political blogosphere.
