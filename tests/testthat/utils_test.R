library(aricode)

rmse <- function(theta, theta_star) { sqrt(sum((theta - theta_star)^2)/sum(theta_star^2)) }

error <- function(beta1, beta2, sort = FALSE) {
  if (sort)
    err <- rmse(sort(beta1), sort(beta2))
  else
    err <- rmse(beta1, beta2)
  err
}

.logistic <- missSBM:::.logistic
.logit <- missSBM:::.logit

# default settings
if (!exists("N_nocov")) N_nocov <- 200
if (!exists("N_cov"))   N_cov <- 40
if (!exists("Q"))       Q <- 3
if (!exists("M"))       M <- 5

## SBM SAMPLER GENERAL PARAMETERS
covarList_directed   <- replicate(M, matrix(rnorm(N_cov * N_cov ,mean = 0, sd = 1), N_cov, N_cov), simplify = FALSE)
covarList_undirected <- lapply(covarList_directed, function(covar) covar + t(covar))
covarParam  <- rnorm(M, 1, 2) * base::sample(c(-1,1), M, replace = TRUE)


## the special case of covariates defined on "nodes"
covarList_node <- replicate(M, rnorm(N_cov ,mean = 0, sd = 1), simplify = FALSE)
covarArray <- missSBM:::getCovarArray(simplify2array(covarList_node), l1_similarity)
covarList_node_similarity <- lapply(1:M, function(m) covarArray[,,m])


## BERNOULLI WITHOUT COVARIATES -------------------------------------------------------------

## UNDIRECTED, NO COVARIATES
sampler_undirected_nocov <- sbm::SimpleSBM$new('bernoulli', N_nocov, FALSE, rep(1/Q, Q), list(mean = diag(.45, Q) + .05 ))
## DIRECTED, NO COVARIATES
sampler_directed_nocov <- sbm::SimpleSBM$new('bernoulli', N_nocov, TRUE, rep(1/Q, Q), list(mean = diag(.45, Q) + matrix(seq(.3, .05, length.out = Q), Q,  Q)))

## BERNOULLI WITH COVARIATES -------------------------------------------------------------

## UNDIRECTED, COVARIATES
sampler_undirected_cov <- sbm::SimpleSBM$new('bernoulli', N_cov, FALSE, rep(1/Q, Q), list(mean = diag(.45, Q) + .05 ), covarParam = covarParam, covarList = covarList_undirected)
sampler_undirected_cov_node <- sbm::SimpleSBM$new('bernoulli', N_cov, FALSE, rep(1/Q, Q), list(mean = diag(.45, Q) + .05 ), covarParam = covarParam, covarList = covarList_node_similarity)

## DIRECTED, COVARIATES
sampler_directed_cov <- sbm::SimpleSBM$new('bernoulli', N_cov, TRUE, rep(1/Q, Q), list(mean = diag(.45, Q) + matrix(seq(.3, .05, length.out = Q), Q,  Q)), covarParam = covarParam, covarList = covarList_directed)

