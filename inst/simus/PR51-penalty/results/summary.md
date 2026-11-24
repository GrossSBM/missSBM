Each cell: all the dyads → the observed dyads in the penalty, mean over 20 replicates.
Commit: 5c8eb2b 



| Network| True Q| Design|   Rate|   n| Observed dyads| True Q selected| Mean selected Q|
|-------:|------:|------:|------:|---:|--------------:|---------------:|---------------:|
|      ER|      1|   dyad|    low| 100|           0.15|     1.00 → 1.00|       1.0 → 1.0|
|      ER|      1|   dyad|    low| 200|           0.15|     1.00 → 1.00|       1.0 → 1.0|
|      ER|      1|   dyad| medium| 100|           0.35|     1.00 → 1.00|       1.0 → 1.0|
|      ER|      1|   dyad| medium| 200|           0.35|     1.00 → 1.00|       1.0 → 1.0|
|      ER|      1|   node|    low| 100|           0.57|     1.00 → 1.00|       1.0 → 1.0|
|      ER|      1|   node|    low| 200|           0.59|     1.00 → 1.00|       1.0 → 1.0|
|      ER|      1|   node| medium| 100|           0.85|     1.00 → 1.00|       1.0 → 1.0|
|      ER|      1|   node| medium| 200|           0.84|     1.00 → 1.00|       1.0 → 1.0|
|   sbm2a|      2|   dyad|    low| 100|           0.15|     0.00 → 0.00|       1.0 → 1.0|
|   sbm2a|      2|   dyad|    low| 200|           0.15|     0.00 → 0.00|       1.0 → 1.0|
|   sbm2a|      2|   dyad| medium| 100|           0.35|     0.00 → 0.00|       1.0 → 1.0|
|   sbm2a|      2|   dyad| medium| 200|           0.35|     0.95 → 0.95|       1.9 → 1.9|
|   sbm2a|      2|   node|    low| 100|           0.58|     0.70 → 0.70|       1.7 → 1.7|
|   sbm2a|      2|   node|    low| 200|           0.57|     1.00 → 1.00|       2.0 → 2.0|
|   sbm2a|      2|   node| medium| 100|           0.84|     1.00 → 1.00|       2.0 → 2.0|
|   sbm2a|      2|   node| medium| 200|           0.84|     1.00 → 1.00|       2.0 → 2.0|
|   sbm2b|      2|   dyad|    low| 100|           0.15|     0.00 → 0.00|       1.0 → 1.0|
|   sbm2b|      2|   dyad|    low| 200|           0.15|     0.10 → 0.10|       1.1 → 1.1|
|   sbm2b|      2|   dyad| medium| 100|           0.35|     0.40 → 0.40|       1.4 → 1.4|
|   sbm2b|      2|   dyad| medium| 200|           0.35|     1.00 → 1.00|       2.0 → 2.0|
|   sbm2b|      2|   node|    low| 100|           0.57|     1.00 → 1.00|       2.0 → 2.0|
|   sbm2b|      2|   node|    low| 200|           0.57|     1.00 → 1.00|       2.0 → 2.0|
|   sbm2b|      2|   node| medium| 100|           0.85|     1.00 → 1.00|       2.0 → 2.0|
|   sbm2b|      2|   node| medium| 200|           0.84|     1.00 → 1.00|       2.0 → 2.0|
|   sbm3a|      3|   dyad|    low| 100|           0.15|     0.00 → 0.00|       1.0 → 1.0|
|   sbm3a|      3|   dyad|    low| 200|           0.15|     0.00 → 0.00|       1.0 → 1.0|
|   sbm3a|      3|   dyad| medium| 100|           0.35|     0.00 → 0.00|       1.0 → 1.0|
|   sbm3a|      3|   dyad| medium| 200|           0.35|     0.75 → 0.80|       2.8 → 2.8|
|   sbm3a|      3|   node|    low| 100|           0.58|     0.10 → 0.10|       1.4 → 1.5|
|   sbm3a|      3|   node|    low| 200|           0.58|     1.00 → 1.00|       3.0 → 3.0|
|   sbm3a|      3|   node| medium| 100|           0.86|     0.80 → 0.80|       2.8 → 2.8|
|   sbm3a|      3|   node| medium| 200|           0.85|     1.00 → 1.00|       3.0 → 3.0|
|   sbm3b|      3|   dyad|    low| 100|           0.15|     0.00 → 0.00|       1.0 → 1.0|
|   sbm3b|      3|   dyad|    low| 200|           0.15|     0.15 → 0.15|       1.4 → 1.5|
|   sbm3b|      3|   dyad| medium| 100|           0.35|     0.10 → 0.20|       1.8 → 1.9|
|   sbm3b|      3|   dyad| medium| 200|           0.35|     1.00 → 1.00|       3.0 → 3.0|
|   sbm3b|      3|   node|    low| 100|           0.58|     0.95 → 0.95|       3.0 → 3.0|
|   sbm3b|      3|   node|    low| 200|           0.58|     1.00 → 1.00|       3.0 → 3.0|
|   sbm3b|      3|   node| medium| 100|           0.84|     1.00 → 1.00|       3.0 → 3.0|
|   sbm3b|      3|   node| medium| 200|           0.84|     1.00 → 1.00|       3.0 → 3.0|


Networks whose selected Q differs between the two penalties: 7 out of 800 



|net   |   n|design |rate   | Qtrue| Q_all| Q_observed|
|:-----|---:|:------|:------|-----:|-----:|----------:|
|sbm3b | 100|dyad   |medium |     3|     2|          3|
|sbm3b | 200|dyad   |low    |     3|     1|          2|
|sbm3b | 200|dyad   |low    |     3|     1|          2|
|sbm3a | 100|node   |low    |     3|     1|          2|
|sbm3b | 100|dyad   |medium |     3|     1|          2|
|sbm3b | 100|dyad   |medium |     3|     2|          3|
|sbm3a | 200|dyad   |medium |     3|     2|          3|
