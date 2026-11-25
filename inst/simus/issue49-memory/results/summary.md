Before (commit c568568) → after (commit da4b5ce). Peaks: increase of the R heap during the phase.



|    n|Design          |partlyObservedNetwork$new()      |clustering()                      |Fit (Q = 3)                       |Total time      |
|----:|:---------------|:--------------------------------|:---------------------------------|:---------------------------------|:---------------|
| 2000|block-dyad      |0.43 s, +190 Mb → 0.16 s, +26 Mb |0.24 s, +205 Mb → 0.14 s, +74 Mb  |1.16 s, +175 Mb → 0.47 s, +67 Mb  |1.8 s → 0.8 s   |
| 4000|block-dyad      |2.78 s, +711 Mb → 0.71 s, +95 Mb |0.66 s, +671 Mb → 0.21 s, +91 Mb  |2.79 s, +637 Mb → 1.18 s, +78 Mb  |6.2 s → 2.1 s   |
| 2000|degree          |0.31 s, +175 Mb → 0.16 s, +26 Mb |0.30 s, +197 Mb → 0.16 s, +90 Mb  |3.44 s, +175 Mb → 3.09 s, +174 Mb |4.0 s → 3.4 s   |
| 4000|degree          |1.23 s, +574 Mb → 0.68 s, +95 Mb |0.64 s, +526 Mb → 0.26 s, +128 Mb |9.98 s, +546 Mb → 9.34 s, +523 Mb |11.8 s → 10.3 s |
| 2000|double-standard |0.44 s, +187 Mb → 0.15 s, +26 Mb |0.26 s, +206 Mb → 0.25 s, +87 Mb  |0.57 s, +204 Mb → 0.76 s, +109 Mb |1.3 s → 1.2 s   |
| 4000|double-standard |2.46 s, +715 Mb → 0.65 s, +95 Mb |0.62 s, +672 Mb → 0.30 s, +117 Mb |1.31 s, +301 Mb → 1.08 s, +110 Mb |4.4 s → 2.0 s   |
| 2000|node            |0.29 s, +221 Mb → 0.18 s, +26 Mb |0.39 s, +244 Mb → 0.28 s, +90 Mb  |0.68 s, +195 Mb → 0.49 s, +3 Mb   |1.4 s → 0.9 s   |
| 4000|node            |1.19 s, +637 Mb → 0.64 s, +95 Mb |0.62 s, +576 Mb → 0.27 s, +127 Mb |1.98 s, +553 Mb → 1.52 s, +4 Mb   |3.8 s → 2.4 s   |
