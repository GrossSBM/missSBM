Each cell: master → PR48-completed → PR48-observed, mean over 20 replicates.
Runs: master at a7cb46d, ICL default (4 min); PR48-completed at e06240d, ICL completed (3 min); PR48-observed at e06240d, ICL observed (3 min).

### Experiment 1: block-node and block-dyad, Q = 1 to 6 explored, true Q = 3

|Design     |Connectivity   |Rate   | Observed dyads|     Q = 3 selected| Mean selected Q|    ARI, selected Q|          ARI, Q = 3| RMSE connectivity, Q = 3|  RMSE sampling, Q = 3| Degenerate fits per collection|
|:----------|:--------------|:------|--------------:|------------------:|---------------:|------------------:|-------------------:|------------------------:|---------------------:|------------------------------:|
|block-dyad |assortative    |low    |           0.16| 0.40 → 0.10 → 1.00| 3.8 → 4.6 → 3.0| 0.61 → 0.81 → 0.99|  0.68 → 0.99 → 0.99|    0.097 → 0.015 → 0.015| 0.062 → 0.052 → 0.049|                0.2 → 0.1 → 0.1|
|block-dyad |assortative    |medium |           0.35| 0.85 → 0.25 → 1.00| 3.2 → 4.2 → 3.0| 0.83 → 0.84 → 0.99|  0.85 → 0.99 → 0.99|    0.059 → 0.011 → 0.011| 0.047 → 0.066 → 0.043|                0.7 → 0.1 → 0.2|
|block-dyad |assortative    |high   |           0.83| 1.00 → 1.00 → 1.00| 3.0 → 3.0 → 3.0| 1.00 → 1.00 → 1.00|  1.00 → 1.00 → 1.00|    0.006 → 0.006 → 0.006| 0.029 → 0.016 → 0.016|                0.7 → 0.1 → 0.1|
|block-dyad |disassortative |low    |           0.16| 0.20 → 0.20 → 0.95| 4.0 → 3.2 → 3.0| 0.00 → 0.30 → 0.97| -0.00 → 0.54 → 0.96|    0.146 → 0.096 → 0.024| 0.097 → 0.091 → 0.036|                0.8 → 0.2 → 0.3|
|block-dyad |disassortative |medium |           0.35| 0.65 → 0.85 → 1.00| 3.0 → 3.2 → 3.0| 0.71 → 0.91 → 0.94|  0.71 → 0.94 → 0.94|    0.050 → 0.020 → 0.019| 0.052 → 0.117 → 0.050|                0.6 → 0.1 → 0.1|
|block-dyad |disassortative |high   |           0.83| 1.00 → 1.00 → 1.00| 3.0 → 3.0 → 3.0| 0.95 → 0.96 → 0.97|  0.95 → 0.96 → 0.97|    0.017 → 0.015 → 0.015| 0.015 → 0.016 → 0.014|                1.4 → 0.1 → 0.1|
|block-node |assortative    |low    |           0.47| 0.10 → 0.10 → 1.00| 4.4 → 5.0 → 3.0| 0.46 → 0.65 → 0.92|  0.38 → 0.62 → 0.92|    0.127 → 0.068 → 0.031| 0.295 → 0.218 → 0.054|                0.2 → 0.0 → 0.0|
|block-node |assortative    |medium |           0.71| 0.00 → 0.05 → 1.00| 4.5 → 4.8 → 3.0| 0.76 → 0.75 → 0.98|  0.78 → 0.93 → 0.98|    0.050 → 0.024 → 0.013| 0.142 → 0.108 → 0.095|                0.2 → 0.0 → 0.0|
|block-node |assortative    |high   |           0.96| 0.80 → 0.50 → 1.00| 3.2 → 3.5 → 3.0| 0.97 → 0.92 → 1.00|  1.00 → 1.00 → 1.00|    0.007 → 0.007 → 0.007| 0.070 → 0.086 → 0.093|                0.8 → 0.0 → 0.1|
|block-node |disassortative |low    |           0.46| 0.00 → 0.40 → 0.15| 3.6 → 3.8 → 2.1| 0.29 → 0.32 → 0.51|  0.19 → 0.28 → 0.49|    0.147 → 0.099 → 0.059| 0.361 → 0.289 → 0.257|                0.7 → 0.0 → 0.0|
|block-node |disassortative |medium |           0.71| 0.35 → 0.10 → 0.95| 4.0 → 4.8 → 3.0| 0.47 → 0.56 → 0.81|  0.44 → 0.58 → 0.81|    0.083 → 0.050 → 0.031| 0.324 → 0.263 → 0.112|                0.8 → 0.0 → 0.0|
|block-node |disassortative |high   |           0.97| 0.95 → 0.65 → 1.00| 3.0 → 3.5 → 3.0| 0.90 → 0.80 → 0.92|  0.91 → 0.85 → 0.92|    0.022 → 0.025 → 0.020| 0.079 → 0.147 → 0.043|                1.4 → 0.0 → 0.0|

### Experiment 2: choice of the design by the ICL (assortative network, medium rate)

|Simulated with  | Observed dyads| True design chosen|Chosen instead (last run)         | Selected Q, true design|   ARI, true design|
|:---------------|--------------:|------------------:|:---------------------------------|-----------------------:|------------------:|
|node            |           0.70| 0.00 → 0.00 → 1.00|block-node: 0, degree: 0          |         3.5 → 3.5 → 3.0| 0.94 → 0.94 → 0.98|
|block-node      |           0.71| 0.85 → 0.80 → 0.95|node: 0, degree: 1                |         4.6 → 4.6 → 3.0| 0.71 → 0.77 → 0.97|
|degree          |           0.67| 0.00 → 0.00 → 0.00|node: 2, block-node: 18           |         4.0 → 3.9 → 3.0| 0.85 → 0.89 → 0.97|
|dyad            |           0.35| 0.95 → 0.90 → 1.00|block-dyad: 0, double-standard: 0 |         5.2 → 5.2 → 2.9| 0.71 → 0.71 → 0.86|
|block-dyad      |           0.35| 0.00 → 0.00 → 1.00|dyad: 0, double-standard: 0       |         3.0 → 3.8 → 3.0| 0.85 → 0.88 → 0.99|
|double-standard |           0.35| 1.00 → 1.00 → 0.60|dyad: 0, block-dyad: 8            |         3.9 → 3.9 → 3.4| 0.84 → 0.84 → 0.93|

Selected Q and ARI of each candidate design, PR48-observed:

|Simulated with  |Fitted with                         |      Selected Q|                ARI|
|:---------------|:-----------------------------------|---------------:|------------------:|
|node            |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.98 / 0.98 / 0.98|
|block-node      |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.97 / 0.97 / 0.97|
|degree          |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.97 / 0.97 / 0.97|
|dyad            |dyad / block-dyad / double-standard | 2.9 / 2.6 / 2.6| 0.86 / 0.75 / 0.67|
|block-dyad      |dyad / block-dyad / double-standard | 3.0 / 3.0 / 3.4| 0.93 / 0.99 / 0.96|
|double-standard |dyad / block-dyad / double-standard | 3.0 / 3.0 / 3.4| 0.98 / 0.99 / 0.93|
