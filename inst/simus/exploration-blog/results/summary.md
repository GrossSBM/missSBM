Runs: master, PR48-completed, PR48-observed. Means over the instances of each group; the ICL gap is measured, for each number of blocks, against the best ICL obtained by the same version over all the settings.
The times are those of runs with 18 jobs in parallel on 20 cores, except `iter2`, run afterwards with 11 jobs: compare them within a run.

### Fully observed network, `node`

|Version        |Setting | Time (s)| ICL gap, mean over Q| ICL gap, max| ICL gap, selected model| Local minima in Q| Collapsed models|Selected Q | ARI with the reference|
|:--------------|:-------|--------:|--------------------:|------------:|-----------------------:|-----------------:|----------------:|:----------|----------------------:|
|master         |none    |      3.2|                406.4|        943.2|                   216.3|               2.0|                1|9          |                   0.67|
|master         |polish  |      3.7|                397.9|        943.2|                   213.6|               2.0|                0|9          |                   0.67|
|master         |default |     22.4|                  1.8|         11.9|                     0.0|               1.0|                0|10         |                   1.00|
|master         |iter3   |     57.7|                  0.0|          0.0|                     0.0|               1.0|                0|10         |                   1.00|
|PR48-completed |none    |      2.6|                406.4|        943.2|                   216.3|               2.0|                1|9          |                   0.67|
|PR48-completed |polish  |      2.8|                397.9|        943.2|                   213.6|               2.0|                0|9          |                   0.67|
|PR48-completed |default |     17.3|                  1.8|         11.9|                     0.0|               1.0|                0|10         |                   1.00|
|PR48-completed |iter3   |     46.1|                  0.0|          0.0|                     0.0|               1.0|                0|10         |                   1.00|
|PR48-observed  |none    |      2.9|                406.4|        943.2|                   216.3|               2.0|                1|9          |                   0.67|
|PR48-observed  |polish  |      3.4|                397.9|        943.2|                   213.6|               2.0|                0|9          |                   0.67|
|PR48-observed  |default |     18.8|                  1.8|         11.9|                     0.0|               1.0|                0|10         |                   1.00|
|PR48-observed  |iter2   |     21.8|                  0.0|          0.0|                     0.0|               1.0|                0|10         |                   1.00|
|PR48-observed  |iter3   |     46.6|                  0.0|          0.0|                     0.0|               1.0|                0|10         |                   1.00|

### 5 block-node samplings, fitted with `node`

|Version        |Setting | Time (s)| ICL gap, mean over Q| ICL gap, max| ICL gap, selected model| Local minima in Q| Collapsed models|Selected Q      | ARI with the reference|
|:--------------|:-------|--------:|--------------------:|------------:|-----------------------:|-----------------:|----------------:|:---------------|----------------------:|
|master         |none    |      2.3|                427.1|       1154.0|                   258.3|               3.8|                3|8, 9, 12, 9, 10 |                   0.56|
|master         |polish  |      3.3|                396.9|       1154.0|                   219.9|               3.6|                0|7, 9, 10, 8, 12 |                   0.60|
|master         |default |     21.2|                 11.2|         98.6|                     6.5|               1.0|                0|9, 9, 9, 9, 10  |                   0.61|
|master         |iter3   |     50.0|                  0.0|          0.0|                     0.0|               1.0|                0|9, 9, 9, 9, 10  |                   0.61|
|PR48-completed |none    |      2.0|                427.1|       1154.0|                   258.3|               3.8|                3|8, 9, 12, 9, 10 |                   0.56|
|PR48-completed |polish  |      2.8|                396.9|       1154.0|                   219.9|               3.6|                0|7, 9, 10, 8, 12 |                   0.60|
|PR48-completed |default |     17.2|                 11.2|         98.6|                     6.5|               1.0|                0|9, 9, 9, 9, 10  |                   0.61|
|PR48-completed |iter3   |     41.9|                  0.0|          0.0|                     0.0|               1.0|                0|9, 9, 9, 9, 10  |                   0.61|
|PR48-observed  |none    |      2.0|                362.1|        940.8|                   174.8|               3.0|                3|8, 7, 9, 9, 7   |                   0.60|
|PR48-observed  |polish  |      2.9|                333.5|        940.8|                   153.8|               2.4|                0|8, 7, 9, 8, 7   |                   0.62|
|PR48-observed  |default |     17.6|                  9.7|        104.0|                     0.5|               1.2|                0|8, 7, 9, 7, 8   |                   0.63|
|PR48-observed  |iter2   |     24.3|                  3.7|         65.7|                     0.0|               1.2|                0|8, 7, 9, 7, 8   |                   0.62|
|PR48-observed  |iter3   |     41.6|                  0.0|          0.0|                     0.0|               1.2|                0|8, 7, 9, 7, 8   |                   0.62|

### 5 block-node samplings, fitted with `block-node`

|Version        |Setting | Time (s)| ICL gap, mean over Q| ICL gap, max| ICL gap, selected model| Local minima in Q| Collapsed models|Selected Q        | ARI with the reference|
|:--------------|:-------|--------:|--------------------:|------------:|-----------------------:|-----------------:|----------------:|:-----------------|----------------------:|
|master         |none    |      4.1|                582.3|       1415.3|                   445.0|               4.2|               60|9, 11, 12, 14, 11 |                   0.00|
|master         |polish  |      5.3|                568.5|       1408.4|                   412.2|               4.0|               59|9, 11, 8, 14, 6   |                   0.12|
|master         |default |     29.7|                 52.3|        753.1|                    16.4|               1.6|                6|8, 7, 7, 11, 8    |                   0.58|
|master         |iter3   |     61.4|                  0.0|          0.0|                     0.0|               1.2|                3|8, 7, 8, 11, 8    |                   0.59|
|PR48-completed |none    |      3.0|                475.9|       1560.0|                   293.2|               3.6|                0|9, 10, 12, 13, 7  |                   0.51|
|PR48-completed |polish  |      5.1|                400.0|       1396.1|                   259.7|               3.8|                0|9, 10, 9, 13, 7   |                   0.53|
|PR48-completed |default |     20.9|                 18.2|        190.0|                    24.3|               2.0|                0|8, 10, 11, 8, 9   |                   0.64|
|PR48-completed |iter3   |     46.3|                  0.0|          0.0|                     0.0|               1.6|                0|8, 11, 11, 8, 9   |                   0.62|
|PR48-observed  |none    |      2.9|                372.9|       1157.4|                   196.2|               3.4|                0|8, 9, 8, 9, 9     |                   0.60|
|PR48-observed  |polish  |      5.7|                312.8|       1157.4|                   152.8|               2.8|                0|8, 9, 9, 7, 9     |                   0.63|
|PR48-observed  |default |     21.3|                  7.6|         95.4|                     3.6|               1.0|                0|8, 7, 9, 8, 8     |                   0.67|
|PR48-observed  |iter2   |     27.7|                  0.1|          9.2|                     0.0|               1.0|                0|8, 7, 9, 8, 9     |                   0.65|
|PR48-observed  |iter3   |     43.7|                  0.0|          0.0|                     0.0|               1.0|                0|8, 7, 9, 8, 9     |                   0.65|
