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
  are node-centered, each entry of `covariates` must be a size-N vector;
  if the covariates are dyad-centered, each entry of `covariates` must
  be N x N matrix.

- similarity:

  an R x R -\> R function to compute similarities between node
  covariates. Default is
  [`l1_similarity`](https://grosssbm.github.io/missSBM/reference/l1_similarity.md),
  that is, -abs(x - y). Only relevant when the covariates are
  node-centered (i.e. `covariates` is a list of size-N vectors).

- intercept:

  An optional intercept term to be added in case of the presence of
  covariates. Default is 0.

## Value

an adjacency matrix with the same dimension as the input, yet with
additional NAs.

## Details

Internal functions use `future_lapply`, so set your plan to
`'multisession'` or `'multicore'` to use several cores/workers.

The sampling designs are either Missing At Random (MAR) or Missing Not
At Random (MNAR), and observe either dyads or nodes (all the dyads of an
observed node). See
[doi:10.1080/01621459.2018.1562934](https://doi.org/10.1080/01621459.2018.1562934)
for a complete description. The parameters of each design:

- Missing At Random (MAR)

  - "dyad": p, the probability that a dyad is observed.

  - "node": p, the probability that a node is observed.

  - "covar-dyad": \\\beta \in R^M\\, a dyad (i, j) being observed with
    probability \\logistic(\beta' x\_{ij})\\, \\x\_{ij}\\ the covariates
    of the dyad.

  - "covar-node": \\\beta \in R^M\\, a node i being observed with
    probability \\logistic(\beta' x_i)\\, \\x_i\\ the covariates of the
    node.

  - "snowball": (number of waves, proportion of the nodes observed in
    the first wave); the nodes of each wave are the neighbors of the
    nodes of the previous one. It is fitted as "node" sampling.

- Missing Not At Random (MNAR)

  - "double-standard": (p0, p1), the probabilities that a dyad is
    observed if it is 0 and if it is 1. (p0, p1) are identified through
    the differences of connectivity between the blocks: with a single
    block, they are not identifiable, and the model is fitted with p0 =
    p1, i.e. as dyad sampling (with a warning).

  - "block-node": \\(p_1, \dots, p_Q)\\, \\p_q\\ the probability that a
    node of block q is observed.

  - "block-dyad": the Q x Q matrix \\(p\_{ql})\\, \\p\_{ql}\\ the
    probability that a dyad between blocks q and l is observed.

  - "degree": (a, b), a node i being observed with probability
    \\logistic(a + b D_i)\\, \\D_i\\ its degree.

The number of blocks is selected by the ICL (Integrated Classification
Likelihood, the smaller the better): -2 times the variational lower
bound of the log-likelihood of the observed data (the observed dyads and
the sampling matrix), at the MAP memberships, plus a penalty, in which
the connectivity counts the observed dyads and the parameters of the
sampling all the dyads (or all the nodes). The ICLs of different designs
fitted to the same data can be compared to choose the design. When their
samplings carry the same information (e.g. "block-node" and "degree" on
an assortative network, where the degree depends on the block), the ICLs
differ by a few units, about the penalty of one parameter (log N), and
the partitions are the same: such a difference is no evidence for one
design. Each number of blocks is fitted from a spectral clustering, then
improved by a polish and by passes of exploration splitting and merging
blocks (see
[`missSBM_param`](https://grosssbm.github.io/missSBM/reference/missSBM_param.md));
one pass is enough on small networks, two are advised for a reference
analysis of a network with many blocks.

## Examples

``` r
## SBM parameters
N <- 300 # number of nodes
Q <- 3   # number of blocks
pi <- rep(1,Q)/Q     # block proportion
theta <- list(mean = diag(.45,Q) + .05 ) # connectivity matrix

## simulate an undirected binary SBM without covariate
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
