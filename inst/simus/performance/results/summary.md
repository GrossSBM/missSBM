v1.1.0 (v1.1.0 installed) → v1.2.0 (v1.2.0 installed). Peaks: increase of the R heap during the phase.



|    n|Design          |partlyObservedNetwork$new()      |clustering()                      |Fit (Q = 3)                         |Total time     |
|----:|:---------------|:--------------------------------|:---------------------------------|:-----------------------------------|:--------------|
| 2000|block-dyad      |0.77 s, +188 Mb → 0.05 s, +27 Mb |0.31 s, +205 Mb → 0.13 s, +74 Mb  |5.78 s, +290 Mb → 0.40 s, +67 Mb    |6.9 s → 0.6 s  |
| 4000|block-dyad      |2.80 s, +711 Mb → 0.21 s, +96 Mb |0.62 s, +671 Mb → 0.18 s, +91 Mb  |8.67 s, +821 Mb → 0.53 s, +77 Mb    |12.1 s → 0.9 s |
| 2000|degree          |0.32 s, +188 Mb → 0.04 s, +27 Mb |0.29 s, +208 Mb → 0.17 s, +90 Mb  |5.29 s, +471 Mb → 1.36 s, +182 Mb   |5.9 s → 1.6 s  |
| 4000|degree          |1.46 s, +716 Mb → 0.12 s, +96 Mb |0.65 s, +660 Mb → 0.21 s, +128 Mb |15.91 s, +1676 Mb → 4.21 s, +452 Mb |18.0 s → 4.5 s |
| 2000|double-standard |0.45 s, +142 Mb → 0.04 s, +27 Mb |0.28 s, +153 Mb → 0.15 s, +87 Mb  |1.80 s, +190 Mb → 0.29 s, +158 Mb   |2.5 s → 0.5 s  |
| 4000|double-standard |2.68 s, +716 Mb → 0.19 s, +96 Mb |0.66 s, +672 Mb → 0.24 s, +117 Mb |5.65 s, +820 Mb → 0.53 s, +165 Mb   |9.0 s → 1.0 s  |
| 2000|node            |0.54 s, +221 Mb → 0.03 s, +27 Mb |0.44 s, +243 Mb → 0.17 s, +91 Mb  |1.08 s, +198 Mb → 0.13 s, +9 Mb     |2.1 s → 0.3 s  |
| 4000|node            |1.19 s, +615 Mb → 0.12 s, +96 Mb |0.66 s, +576 Mb → 0.21 s, +128 Mb |2.90 s, +556 Mb → 0.39 s, +10 Mb    |4.7 s → 0.7 s  |
