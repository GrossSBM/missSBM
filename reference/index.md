# Package index

## Top-level fitting functions

Main functions for estimating and sampling from an SBM with missing data

- [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)
  : Estimation of simple SBMs with missing data
- [`observeNetwork()`](https://grosssbm.github.io/missSBM/reference/observeNetwork.md)
  : Observe a network partially according to a given sampling design
- [`missSBM_param()`](https://grosssbm.github.io/missSBM/reference/missSBM_param.md)
  : Control of a missSBM fit

## Main classes of objects

Description of objects missSBM_fit and missSBM_collection. The class
missSBM_fit is the more central class of object, embedding fits for both
the SBM and the sampling model. The class missSBM_collection defines
objects for storing a collection of missSBM_fit, resulting from the the
top-level function estimateMissSBM().

- [`missSBM_fit`](https://grosssbm.github.io/missSBM/reference/missSBM_fit.md)
  : An R6 class to represent an SBM fit with missing data

- [`coef(`*`<missSBM_fit>`*`)`](https://grosssbm.github.io/missSBM/reference/coef.missSBM_fit.md)
  : Extract model coefficients

- [`fitted(`*`<missSBM_fit>`*`)`](https://grosssbm.github.io/missSBM/reference/fitted.missSBM_fit.md)
  :

  Extract model fitted values from object `missSBM_fit`, return by
  [`estimateMissSBM()`](https://grosssbm.github.io/missSBM/reference/estimateMissSBM.md)

- [`predict(`*`<missSBM_fit>`*`)`](https://grosssbm.github.io/missSBM/reference/predicted.missSBM_fit.md)
  :

  Prediction of a `missSBM_fit` (i.e. network with imputed missing
  dyads)

- [`plot(`*`<missSBM_fit>`*`)`](https://grosssbm.github.io/missSBM/reference/plot.missSBM_fit.md)
  :

  Visualization for an object `missSBM_fit`

- [`missSBM_collection`](https://grosssbm.github.io/missSBM/reference/missSBM_collection.md)
  : An R6 class to represent a collection of SBM fits with missing data

## Data sets

- [`war`](https://grosssbm.github.io/missSBM/reference/war.md) : War
  data set
- [`frenchblog2007`](https://grosssbm.github.io/missSBM/reference/frenchblog2007.md)
  : Political Blogosphere network prior to 2007 French presidential
  election
- [`er_network`](https://grosssbm.github.io/missSBM/reference/er_network.md)
  : ER ego centered network
