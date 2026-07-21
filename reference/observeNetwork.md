# Observe a network partially according to a given sampling design

This function draws observations in an adjacency matrix according to a
given network sampling design.

## Usage

``` r
observeNetwork(
  adjacencyMatrix,
  sampling,
  parameters,
  clusters = NULL,
  covariates = list(),
  similarity = l1_similarity,
  intercept = 0
)
```

## Arguments

- adjacencyMatrix:

  The N x N adjacency matrix of the network to sample. The diagonal is
  expected to be NA (no self-loops); any other pre-existing NA entry is
  treated as an absent edge (coded 0) before sampling is applied on top
  of it, with a warning.

- sampling:

  The sampling design used to observe the adjacency matrix, see details.

- parameters:

  The sampling parameters (adapted to each sampling, see details).

- clusters:

  An optional clustering membership vector of the nodes. Only necessary
  for block samplings.

- covariates:

  An optional list with M entries (the M covariates). If the covariates
  are node-centered, each entry of `covariates`. must be a size-N
  vector; if the covariates are dyad-centered, each entry of
  `covariates` must be N x N matrix.

- similarity:

  An optional function to compute similarities between node covariates.
  Default is
  [`l1_similarity`](https://grosssbm.github.io/missSBM/reference/l1_similarity.md),
  that is, -abs(x-y). Only relevant when the covariates are
  node-centered.

- intercept:

  An optional intercept term to be added in case of the presence of
  covariates. Default is 0.

## Value

an adjacency matrix with the same dimension as the input, yet with
additional NAs.

## Details

Internal functions use `future_lapply`, so set your plan to
`'multisession'` or `'multicore'` to use several cores/workers.

The different sampling designs are split into two families in which we
find dyad-centered and node-centered samplings. See
[doi:10.1080/01621459.2018.1562934](https://doi.org/10.1080/01621459.2018.1562934)
for a complete description.

- Missing at Random (MAR)

  - dyad parameter = p = Prob(Dyad(i,j) is observed)

  - node parameter = p = Prob(Node i is observed)

  - covar-dyad": parameter = beta in R^M, such that Prob(Dyad (i,j) is
    observed) = logistic(parameter' covarArray (i,j, .))

  - covar-node": parameter = nu in R^M such that Prob(Node i is
    observed) = logistic(parameter' covarMatrix (i,)

  - snowball": parameter = number of waves with Prob(Node i is observed
    in the 1st wave)

- Missing Not At Random (MNAR)

  - double-standard parameter = (p0,p1) with p0 = Prob(Dyad (i,j) is
    observed \| the dyad is equal to 0), p1 = Prob(Dyad (i,j) is
    observed \| the dyad is equal to 1)

  - block-node parameter = c(p(1),...,p(Q)) and p(q) = Prob(Node i is
    observed \| node i is in cluster q)

  - block-dyad parameter = c(p(1,1),...,p(Q,Q)) and p(q,l) = Prob(Edge
    (i,j) is observed \| node i is in cluster q and node j is in cluster
    l)

## Examples

``` r
## SBM parameters
N <- 300 # number of nodes
Q <- 3   # number of clusters
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix

## simulate an unidrected binary SBM without covariate
sbm <- sbm::sampleSimpleSBM(N, pi, theta)

## Sample network data

# some sampling design and their associated parameters
sampling_parameters <- list(
   "dyad" = .3,
   "node" = .3,
   "double-standard" = c(0.4, 0.8),
   "block-node" = c(.3, .8, .5),
   "block-dyad" = theta$mean,
   "degree" = c(.01, .01),
   "snowball" = c(2,.1)
 )

observed_networks <- list()

for (sampling in names(sampling_parameters)) {
  observed_networks[[sampling]] <-
     missSBM::observeNetwork(
       adjacencyMatrix = sbm$networkData,
       sampling        = sampling,
       parameters      = sampling_parameters[[sampling]],
       clusters        = sbm$memberships
     )
}
```
