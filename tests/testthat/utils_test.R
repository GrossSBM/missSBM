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

# N_nocov <- 200
# N_cov   <- 40
# Q <- 3
# M <- 5

## SBM SAMPLER GENERAL PARAMETERS
covarList_directed   <- replicate(M, matrix(rnorm(N_cov * N_cov ,mean = 0, sd = 1), N_cov, N_cov), simplify = FALSE)
covarList_undirected <- lapply(covarList_directed, function(covar) covar + t(covar))
covarParam  <- rnorm(M, 1, 2) * base::sample(c(-1,1), M, replace = TRUE)

## BERNOULLI WITHOUT COVARIATES -------------------------------------------------------------

## UNDIRECTED, NO COVARIATES
sampler_undirected_nocov <- sbm::SimpleSBM$new('bernoulli', N_nocov, FALSE, rep(1/Q, Q), list(mean = diag(.45, Q) + .05 ))
## DIRECTED, NO COVARIATES
sampler_directed_nocov <- sbm::SimpleSBM$new('bernoulli', N_nocov, TRUE, rep(1/Q, Q), list(mean = matrix(seq(.9, .1, length.out = Q*Q), Q,  Q)))

## BERNOULLI WITH COVARIATES -------------------------------------------------------------

## UNDIRECTED, COVARIATES
sampler_undirected_cov <- sbm::SimpleSBM$new('bernoulli', N_cov, FALSE, rep(1/Q, Q), list(mean = diag(.45, Q) + .05 ), covarParam = covarParam, covarList = covarList_undirected)

## DIRECTED, COVARIATES
sampler_directed_cov <- sbm::SimpleSBM$new('bernoulli', N_cov, TRUE, rep(1/Q, Q), list(mean = matrix(seq(.9, .1, length.out = Q*Q), Q,  Q)), covarParam = covarParam, covarList = covarList_directed)

