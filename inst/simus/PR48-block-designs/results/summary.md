Each cell: master → PR48, mean over 20 replicates. Runs: master at a7cb46d (3 min), PR48 at 4802f93 (2 min).
The last columns: PR48 only, models selected with the ICL built on the observed data (see README).

### Experiment 1: block-node and block-dyad, Q = 1 to 6 explored, true Q = 3

|Design     |Connectivity   |Rate   | Observed dyads| Q = 3 selected| Mean selected Q| ARI, selected Q|   ARI, Q = 3| RMSE connectivity, Q = 3| RMSE sampling, Q = 3| Degenerate fits per collection| Q = 3 selected, observed-data ICL| ARI, observed-data ICL|
|:----------|:--------------|:------|--------------:|--------------:|---------------:|---------------:|------------:|------------------------:|--------------------:|------------------------------:|---------------------------------:|----------------------:|
|block-dyad |assortative    |low    |           0.16|    0.40 → 0.10|       3.8 → 4.6|     0.61 → 0.81|  0.68 → 0.99|            0.097 → 0.015|        0.062 → 0.052|                      0.2 → 0.1|                              1.00|                   0.99|
|block-dyad |assortative    |medium |           0.35|    0.85 → 0.25|       3.2 → 4.2|     0.83 → 0.84|  0.85 → 0.99|            0.059 → 0.011|        0.047 → 0.066|                      0.7 → 0.1|                              1.00|                   0.99|
|block-dyad |assortative    |high   |           0.83|    1.00 → 1.00|       3.0 → 3.0|     1.00 → 1.00|  1.00 → 1.00|            0.006 → 0.006|        0.029 → 0.016|                      0.7 → 0.1|                              1.00|                   1.00|
|block-dyad |disassortative |low    |           0.16|    0.20 → 0.20|       4.0 → 3.2|     0.00 → 0.30| -0.00 → 0.54|            0.146 → 0.096|        0.097 → 0.091|                      0.8 → 0.2|                              0.20|                   0.77|
|block-dyad |disassortative |medium |           0.35|    0.65 → 0.85|       3.0 → 3.2|     0.71 → 0.91|  0.71 → 0.94|            0.050 → 0.020|        0.052 → 0.117|                      0.6 → 0.1|                              1.00|                   0.94|
|block-dyad |disassortative |high   |           0.83|    1.00 → 1.00|       3.0 → 3.0|     0.95 → 0.96|  0.95 → 0.96|            0.017 → 0.015|        0.015 → 0.016|                      1.4 → 0.1|                              1.00|                   0.96|
|block-node |assortative    |low    |           0.47|    0.10 → 0.10|       4.4 → 5.0|     0.46 → 0.65|  0.38 → 0.62|            0.127 → 0.068|        0.295 → 0.218|                      0.2 → 0.0|                              0.30|                   0.74|
|block-node |assortative    |medium |           0.71|    0.00 → 0.05|       4.5 → 4.8|     0.76 → 0.75|  0.78 → 0.93|            0.050 → 0.024|        0.142 → 0.108|                      0.2 → 0.0|                              0.85|                   0.96|
|block-node |assortative    |high   |           0.96|    0.80 → 0.50|       3.2 → 3.5|     0.97 → 0.92|  1.00 → 1.00|            0.007 → 0.007|        0.070 → 0.086|                      0.8 → 0.0|                              1.00|                   1.00|
|block-node |disassortative |low    |           0.46|    0.00 → 0.40|       3.6 → 3.8|     0.29 → 0.32|  0.19 → 0.28|            0.147 → 0.099|        0.361 → 0.289|                      0.7 → 0.0|                              0.30|                   0.45|
|block-node |disassortative |medium |           0.71|    0.35 → 0.10|       4.0 → 4.8|     0.47 → 0.56|  0.44 → 0.58|            0.083 → 0.050|        0.324 → 0.263|                      0.8 → 0.0|                              0.50|                   0.58|
|block-node |disassortative |high   |           0.97|    0.95 → 0.65|       3.0 → 3.5|     0.90 → 0.80|  0.91 → 0.85|            0.022 → 0.025|        0.079 → 0.147|                      1.4 → 0.0|                              0.85|                   0.90|

### Experiment 2: choice of the design by the ICL (assortative network, medium rate)

|Simulated with | Observed dyads| True design chosen| Selected Q, MAR design| Selected Q, block design| ARI, MAR design| ARI, block design| True design chosen, observed-data ICL| Selected Q (MAR / block), observed-data ICL|
|:--------------|--------------:|------------------:|----------------------:|------------------------:|---------------:|-----------------:|-------------------------------------:|-------------------------------------------:|
|node           |           0.70|        0.00 → 0.00|              3.5 → 3.5|                4.8 → 5.0|     0.94 → 0.94|       0.69 → 0.66|                                  1.00|                                   3.0 / 3.6|
|block-node     |           0.71|        0.90 → 0.80|              3.9 → 3.9|                4.6 → 4.6|     0.90 → 0.90|       0.71 → 0.77|                                  0.95|                                   3.0 / 3.0|
|dyad           |           0.35|        1.00 → 1.00|              5.5 → 5.5|                3.5 → 4.5|     0.67 → 0.67|       0.69 → 0.72|                                  1.00|                                   3.0 / 2.8|
|block-dyad     |           0.35|        0.95 → 1.00|              5.0 → 5.0|                3.2 → 3.8|     0.75 → 0.75|       0.80 → 0.89|                                  1.00|                                   3.0 / 3.0|
