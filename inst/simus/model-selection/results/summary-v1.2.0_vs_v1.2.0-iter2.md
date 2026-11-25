Each cell: v1.2.0 → v1.2.0-iter2, mean over 20 replicates.
Runs: v1.2.0, installed (2 min on 16 cores); v1.2.0-iter2, installed (4 min on 16 cores).

### Experiment 1: block-node and block-dyad, Q = 1 to 6 explored, true Q = 3

|Design     |Connectivity   |Rate   | Observed dyads| Q = 3 selected| Mean selected Q| ARI, selected Q|  ARI, Q = 3| RMSE connectivity, Q = 3| RMSE sampling, Q = 3| Degenerate fits per collection|
|:----------|:--------------|:------|--------------:|--------------:|---------------:|---------------:|-----------:|------------------------:|--------------------:|------------------------------:|
|block-dyad |assortative    |low    |           0.16|    0.95 → 0.95|       3.0 → 3.0|     0.97 → 0.97| 0.97 → 0.97|            0.018 → 0.018|        0.012 → 0.012|                      0.2 → 0.1|
|block-dyad |assortative    |medium |           0.35|    1.00 → 1.00|       3.0 → 3.0|     0.99 → 0.99| 0.99 → 0.99|            0.011 → 0.011|        0.012 → 0.012|                      0.1 → 0.1|
|block-dyad |assortative    |high   |           0.83|    1.00 → 1.00|       3.0 → 3.0|     1.00 → 1.00| 1.00 → 1.00|            0.006 → 0.006|        0.013 → 0.013|                      0.1 → 0.1|
|block-dyad |disassortative |low    |           0.16|    0.85 → 0.85|       3.1 → 3.1|     0.94 → 0.94| 0.89 → 0.89|            0.032 → 0.032|        0.021 → 0.021|                      0.3 → 0.3|
|block-dyad |disassortative |medium |           0.35|    1.00 → 1.00|       3.0 → 3.0|     0.94 → 0.94| 0.94 → 0.94|            0.019 → 0.019|        0.021 → 0.021|                      0.1 → 0.1|
|block-dyad |disassortative |high   |           0.83|    1.00 → 1.00|       3.0 → 3.0|     0.97 → 0.97| 0.97 → 0.97|            0.014 → 0.014|        0.014 → 0.014|                      0.1 → 0.1|
|block-node |assortative    |low    |           0.47|    1.00 → 1.00|       3.0 → 3.0|     0.92 → 0.92| 0.92 → 0.92|            0.031 → 0.031|        0.054 → 0.054|                      0.0 → 0.0|
|block-node |assortative    |medium |           0.71|    1.00 → 1.00|       3.0 → 3.0|     0.99 → 0.99| 0.99 → 0.99|            0.013 → 0.013|        0.050 → 0.050|                      0.0 → 0.0|
|block-node |assortative    |high   |           0.96|    1.00 → 1.00|       3.0 → 3.0|     1.00 → 1.00| 1.00 → 1.00|            0.007 → 0.007|        0.055 → 0.055|                      0.0 → 0.0|
|block-node |disassortative |low    |           0.46|    0.25 → 0.25|       2.2 → 2.2|     0.53 → 0.53| 0.57 → 0.57|            0.054 → 0.054|        0.150 → 0.150|                      0.0 → 0.0|
|block-node |disassortative |medium |           0.71|    0.95 → 0.95|       3.0 → 3.0|     0.81 → 0.81| 0.82 → 0.82|            0.030 → 0.030|        0.089 → 0.089|                      0.0 → 0.0|
|block-node |disassortative |high   |           0.97|    1.00 → 1.00|       3.0 → 3.0|     0.92 → 0.92| 0.92 → 0.92|            0.020 → 0.020|        0.043 → 0.043|                      0.0 → 0.0|

### Experiment 2: choice of the design by the ICL (assortative network, medium rate)

|Simulated with  | Observed dyads| True design chosen|Chosen instead (last run)         | Selected Q, true design| ARI, true design|
|:---------------|--------------:|------------------:|:---------------------------------|-----------------------:|----------------:|
|node            |           0.70|        0.95 → 0.95|block-node: 0, degree: 1          |               3.0 → 3.0|      0.98 → 0.98|
|block-node      |           0.71|        0.60 → 0.60|node: 0, degree: 8                |               3.0 → 3.0|      0.97 → 0.97|
|degree          |           0.67|        0.90 → 0.90|node: 1, block-node: 1            |               3.0 → 3.0|      0.97 → 0.97|
|dyad            |           0.35|        1.00 → 1.00|block-dyad: 0, double-standard: 0 |               3.0 → 3.0|      0.90 → 0.90|
|block-dyad      |           0.35|        1.00 → 1.00|dyad: 0, double-standard: 0       |               3.0 → 3.0|      0.99 → 0.99|
|double-standard |           0.35|        1.00 → 1.00|dyad: 0, block-dyad: 0            |               3.0 → 3.0|      0.99 → 0.99|

Selected Q and ARI of each candidate design, v1.2.0-iter2:

|Simulated with  |Fitted with                         |      Selected Q|                ARI|
|:---------------|:-----------------------------------|---------------:|------------------:|
|node            |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.98 / 0.98 / 0.98|
|block-node      |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.97 / 0.97 / 0.97|
|degree          |node / block-node / degree          | 3.0 / 3.0 / 3.0| 0.97 / 0.97 / 0.97|
|dyad            |dyad / block-dyad / double-standard | 3.0 / 2.6 / 3.0| 0.90 / 0.77 / 0.91|
|block-dyad      |dyad / block-dyad / double-standard | 3.0 / 3.0 / 3.0| 0.93 / 0.99 / 0.98|
|double-standard |dyad / block-dyad / double-standard | 3.0 / 3.0 / 3.0| 0.98 / 0.99 / 0.99|
