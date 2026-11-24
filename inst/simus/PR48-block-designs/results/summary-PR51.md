Each cell: master → PR52-observed → PR51-observed, mean over 20 replicates.
Runs: master at a7cb46d, ICL default (4 min); PR52-observed at cbe717a, ICL observed (4 min); PR51-observed at 522d354, ICL default (4 min).

### Experiment 1: block-node and block-dyad, Q = 1 to 6 explored, true Q = 3

|Design     |Connectivity   |Rate   | Observed dyads|     Q = 3 selected| Mean selected Q|    ARI, selected Q|          ARI, Q = 3| RMSE connectivity, Q = 3|  RMSE sampling, Q = 3| Degenerate fits per collection|
|:----------|:--------------|:------|--------------:|------------------:|---------------:|------------------:|-------------------:|------------------------:|---------------------:|------------------------------:|
|block-dyad |assortative    |low    |           0.16| 0.40 → 0.95 → 0.95| 3.8 → 3.0 → 3.0| 0.61 → 0.97 → 0.97|  0.68 → 0.97 → 0.97|    0.097 → 0.018 → 0.018| 0.062 → 0.012 → 0.012|                0.2 → 0.3 → 0.3|
|block-dyad |assortative    |medium |           0.35| 0.85 → 1.00 → 1.00| 3.2 → 3.0 → 3.0| 0.83 → 0.99 → 0.99|  0.85 → 0.99 → 0.99|    0.059 → 0.011 → 0.011| 0.047 → 0.012 → 0.012|                0.7 → 0.2 → 0.2|
|block-dyad |assortative    |high   |           0.83| 1.00 → 1.00 → 1.00| 3.0 → 3.0 → 3.0| 1.00 → 1.00 → 1.00|  1.00 → 1.00 → 1.00|    0.006 → 0.006 → 0.006| 0.029 → 0.013 → 0.013|                0.7 → 0.1 → 0.1|
|block-dyad |disassortative |low    |           0.16| 0.20 → 0.85 → 0.85| 4.0 → 3.1 → 3.1| 0.00 → 0.94 → 0.94| -0.00 → 0.89 → 0.89|    0.146 → 0.032 → 0.032| 0.097 → 0.021 → 0.021|                0.8 → 0.3 → 0.3|
|block-dyad |disassortative |medium |           0.35| 0.65 → 1.00 → 1.00| 3.0 → 3.0 → 3.0| 0.71 → 0.94 → 0.94|  0.71 → 0.94 → 0.94|    0.050 → 0.019 → 0.019| 0.052 → 0.021 → 0.021|                0.6 → 0.1 → 0.1|
|block-dyad |disassortative |high   |           0.83| 1.00 → 1.00 → 1.00| 3.0 → 3.0 → 3.0| 0.95 → 0.96 → 0.96|  0.95 → 0.96 → 0.96|    0.017 → 0.015 → 0.015| 0.015 → 0.015 → 0.015|                1.4 → 0.1 → 0.1|
|block-node |assortative    |low    |           0.47| 0.10 → 1.00 → 1.00| 4.4 → 3.0 → 3.0| 0.46 → 0.92 → 0.92|  0.38 → 0.92 → 0.92|    0.127 → 0.031 → 0.031| 0.295 → 0.054 → 0.054|                0.2 → 0.0 → 0.0|
|block-node |assortative    |medium |           0.71| 0.00 → 1.00 → 1.00| 4.5 → 3.0 → 3.0| 0.76 → 0.98 → 0.98|  0.78 → 0.98 → 0.98|    0.050 → 0.013 → 0.013| 0.142 → 0.050 → 0.050|                0.2 → 0.0 → 0.0|
|block-node |assortative    |high   |           0.96| 0.80 → 1.00 → 1.00| 3.2 → 3.0 → 3.0| 0.97 → 1.00 → 1.00|  1.00 → 1.00 → 1.00|    0.007 → 0.007 → 0.007| 0.070 → 0.055 → 0.055|                0.8 → 0.0 → 0.0|
|block-node |disassortative |low    |           0.46| 0.00 → 0.15 → 0.15| 3.6 → 2.1 → 2.1| 0.29 → 0.49 → 0.49|  0.19 → 0.47 → 0.47|    0.147 → 0.060 → 0.060| 0.361 → 0.213 → 0.213|                0.7 → 0.0 → 0.0|
|block-node |disassortative |medium |           0.71| 0.35 → 0.95 → 0.95| 4.0 → 3.0 → 3.0| 0.47 → 0.81 → 0.81|  0.44 → 0.80 → 0.80|    0.083 → 0.031 → 0.031| 0.324 → 0.101 → 0.101|                0.8 → 0.0 → 0.0|
|block-node |disassortative |high   |           0.97| 0.95 → 1.00 → 1.00| 3.0 → 3.0 → 3.0| 0.90 → 0.92 → 0.92|  0.91 → 0.92 → 0.92|    0.022 → 0.020 → 0.020| 0.079 → 0.043 → 0.043|                1.4 → 0.0 → 0.0|

### Experiment 2: choice of the design by the ICL (assortative network, medium rate)

|Simulated with  | Observed dyads| True design chosen|Chosen instead (last run)         | Selected Q, true design|   ARI, true design|
|:---------------|--------------:|------------------:|:---------------------------------|-----------------------:|------------------:|
|node            |           0.70| 0.00 → 1.00 → 1.00|block-node: 0, degree: 0          |         3.5 → 3.0 → 3.0| 0.94 → 0.98 → 0.98|
|block-node      |           0.71| 0.85 → 0.95 → 0.95|node: 0, degree: 1                |         4.6 → 3.0 → 3.0| 0.71 → 0.97 → 0.97|
|degree          |           0.67| 0.00 → 0.00 → 0.00|node: 2, block-node: 18           |         4.0 → 3.0 → 3.0| 0.85 → 0.97 → 0.97|
|dyad            |           0.35| 0.95 → 1.00 → 1.00|block-dyad: 0, double-standard: 0 |         5.2 → 2.9 → 3.0| 0.71 → 0.86 → 0.89|
|block-dyad      |           0.35| 0.00 → 1.00 → 1.00|dyad: 0, double-standard: 0       |         3.0 → 3.0 → 3.0| 0.85 → 0.99 → 0.99|
|double-standard |           0.35| 1.00 → 1.00 → 1.00|dyad: 0, block-dyad: 0            |         3.9 → 3.0 → 3.0| 0.84 → 0.99 → 0.99|

Selected Q and ARI of each candidate design, PR51-observed:

|Simulated with  |Fitted with                         |      Selected Q|                ARI|
|:---------------|:-----------------------------------|---------------:|------------------:|
|node            |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.98 / 0.98 / 0.98|
|block-node      |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.97 / 0.97 / 0.97|
|degree          |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.97 / 0.97 / 0.97|
|dyad            |dyad / block-dyad / double-standard | 3.0 / 2.6 / 3.0| 0.89 / 0.75 / 0.89|
|block-dyad      |dyad / block-dyad / double-standard | 3.0 / 3.0 / 3.0| 0.93 / 0.99 / 0.98|
|double-standard |dyad / block-dyad / double-standard | 3.0 / 3.0 / 3.0| 0.98 / 0.99 / 0.99|
