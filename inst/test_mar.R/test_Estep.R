library(sbm)
N_cov <- 100
Q <- 2
M <- 10
source("tests/testthat/utils_test.R")

sampler_undirected_cov$rNetwork(store = TRUE)
net <- missSBM:::partlyObservedNetwork$new(sampler_undirected_cov$networkData, covariates = covarList_undirected)
cls <- net$clustering(1:(2*Q))
sbm <- missSBM:::SimpleSBM_fit_withCov$new(net, clusterInit = cls[[2]], covarList = covarList_undirected)

Gamma <- .logit(sbm$connectParam$mean)
pi <- sbm$blockProp
Tau <- sbm$probMemberships

spmat = missSBM:::E_step_sparse_bernoulli_covariates(net$networkData, net$samplingMatrix, sampler_undirected_cov$covarEffect, Tau, Gamma, pi, TRUE)
dense = missSBM:::E_step_covariates(sampler_undirected_cov$networkData, sampler_undirected_cov$covarEffect, Gamma, Tau, pi)

plot(spmat, dense)

cl_sp <- apply(spmat, 1, which.max)
cl_ds <- apply(dense, 1, which.max)
aricode::ARI(cls[[Q]], sampler_undirected_cov$memberships)
aricode::ARI(cl_ds, cl_sp)
aricode::ARI(cl_sp, sampler_undirected_cov$memberships)
aricode::ARI(cl_ds, sampler_undirected_cov$memberships)

microbenchmark::microbenchmark(
spmat = missSBM:::E_step_sparse_bernoulli_covariates(net$networkData, net$samplingMatrix, sampler_undirected_cov$covarEffect, Tau, Gamma, pi, TRUE),
dense = missSBM:::E_step_covariates(sampler_undirected_cov$networkData, sampler_undirected_cov$covarEffect, Gamma, Tau, pi)
) -> res_Estep

ggplot2::autoplot(res_Estep)

