library(sbm)
library(missSBM)

## Common parameters
nbNodes  <- 500
nbBlocks <- 2
blockProp <- c(.5, .5) # group proportions
covarParam <- c(-2,2)
dimLabels <- list(row = "rowLabel", col = "colLabel")
covar1 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covar2 <- matrix(rnorm(nbNodes**2), nbNodes, nbNodes)
covarList_directed <- list(covar1 = covar1, covar2 = covar2)

covar1 <- covar1 + t(covar1)
covar2 <- covar2 + t(covar2)
covarList <- list(covar1 = covar1, covar2 = covar2)

means <- diag(.4, 2) + 0.05
connectParam <- list(mean = means)

## Basic construction - check for wrong specifications
mySampler <- sbm::SimpleSBM$new('bernoulli', nbNodes, FALSE, blockProp, connectParam, covarParam = covarParam[1], covarList = covarList[1])
mySampler$rMemberships(store = TRUE)
mySampler$rEdges(store = TRUE)
adjacencyMatrix <- missSBM::observeNetwork(mySampler$networkData, "dyad", 1)
A <- adjacencyMatrix

Y <- adjacencyMatrix
Y[is.na(Y)] <- 0
Z <- mySampler$indMemberships
Tau <- matrix(0.5, nrow(Z), ncol(Z))
M <- mySampler$covarEffect
Gamma <- missSBM:::.logit(mySampler$connectParam$mean)
pi <- mySampler$blockProp

# Storing data
## where are my observations?
diag(A) <- NA
obs <- which(!is.na(A), arr.ind = TRUE)
R_sp  <- Matrix::sparseMatrix(obs[,1], obs[,2],x = 1, dims = dim(A))
obs   <- which(!is.na(A) & upper.tri(A), arr.ind = TRUE)
R_sp2 <- Matrix::sparseMatrix(obs[,1], obs[,2],x = 1, dims = dim(A))

## where are my non-zero entries?
nzero <- which(!is.na(A) & A != 0 , arr.ind = TRUE)
Y_sp   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(A))
nzero <- which(!is.na(A) & A != 0  & upper.tri(A), arr.ind = TRUE)
Y_sp2   <- Matrix::sparseMatrix(nzero[,1], nzero[,2], x = 1, dims = dim(A))

microbenchmark::microbenchmark(
 spmat = missSBM:::vLL_complete_sparse_bernoulli_covariates(Y_sp2, R_sp2, M, Tau, Gamma, pi),
 dense = missSBM:::vExpec_covariates(Y, M, Gamma, Tau, pi)
) -> res_vLL


microbenchmark::microbenchmark(
 spmat = missSBM:::E_step_sparse_bernoulli_covariates(Y_sp2, R_sp2, M, Tau, Gamma, pi),
 dense = missSBM:::E_step_covariates(Y, M, Gamma, Tau, pi)
) -> res_Estep

# Z1 <- missSBM:::E_step_sparse_bernoulli_undirected_covariates(Y_sp, R_sp, M, Z, Gamma, pi)
# Z2 <- missSBM:::E_step_sparse_bernoulli_covariates(Y_sp2, R_sp2, M, Z, Gamma, pi)

param <- list(Gamma = Gamma, beta = matrix(covarParam[1], ncol=1))

microbenchmark::microbenchmark(
sparse = sparse <- missSBM:::M_step_sparse_bernoulli_covariates(
        param,
        Y_sp2,
        R_sp2,
        mySampler$covarArray,
        Tau,
        configuration = list(algorithm="CCSAQ", maxeval = 50, xtol_rel = 1e-4)
        ),
dense  = nloptr::nloptr(
                # starting parameters
                unlist(param),
                # objective function + gradient
                missSBM:::Mstep_covariates_undirected,
                # optimizer parameters
                opts = list("algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-4),
                # additional argument for objective/gradient function
                Y = Y, cov = mySampler$covarArray, Tau = Tau,
        ), times = 10) -> res_Mstep_covariates

ggplot2::autoplot(res_vLL)
ggplot2::autoplot(res_Estep)
ggplot2::autoplot(res_Mstep_covariates)

beta_new <- sparse$beta
gamma_new <- matrix(.logit(sparse$theta$mean), 2, 2)
beta_old  <- dense$solution[-(1:4)]
gamma_old <- matrix(dense$solution[1:4], 2, 2)
