## Issue #55: the three estimators of the degree sampling (see proto.R) on networks simulated
## under the degree sampling, with an effect of the degree from none to strong.
##
## Usage, from the root of the package:
##   OPENBLAS_NUM_THREADS=1 Rscript inst/simus/issue55-degree/simulate.R [cores] [n]
args  <- commandArgs(TRUE)
cores <- if (length(args) > 0) as.integer(args[1]) else 16
n     <- if (length(args) > 1) as.integer(args[2]) else 200
suppressMessages(devtools::load_all(quiet = TRUE))
source("inst/simus/issue55-degree/proto.R")
future::plan("sequential")
library(parallel)

## the network of the campaign of #48 (assortative), Q = 3
Q     <- 3
alpha <- c(.25, .35, .40)
pi    <- matrix(.05, Q, Q); diag(pi) <- c(.6, .3, .15)
## (a, b): about 43% of observed nodes in all the settings, the effect of the degree increasing;
## b is given for n = 200 and scaled by 200 / n, the degrees growing with n
psis  <- lapply(list(none = c(-.2, 0), weak = c(-1.4, .04), medium = c(-2.4, .08), strong = c(-4, .15)),
                function(p) c(p[1], p[2] * 200 / n))
grid  <- expand.grid(effect = names(psis), rep = 1:20, stringsAsFactors = FALSE)

## RMSE of the connectivity at the level of the dyads (invariant by relabelling)
rmse_pairs <- function(cl, theta, cl0) {
  up <- upper.tri(matrix(0, length(cl), length(cl)))
  sqrt(mean((theta[cl, cl][up] - pi[cl0, cl0][up])^2))
}

run <- function(k) {
  g <- grid[k, ]; set.seed(k)
  s0  <- sbm::sampleSimpleSBM(n, alpha, list(mean = pi))
  psi <- psis[[g$effect]]
  A   <- observeNetwork(s0$networkData, "degree", psi)
  pon <- partlyObservedNetwork$new(A)
  cl  <- pon$clustering(Q)[[1]]

  t_A <- system.time(stA <- fit_A(pon, cl))[3]
  t_B <- system.time(stB <- fit_B(pon, cl))[3]
  t_C <- system.time(stC <- fit_C(stA))[3] + t_A
  ## the other designs of the family, with the package (at the true Q, no exploration)
  ctrl <- missSBM_param(trace = FALSE, iterates = 0)
  ICL_node       <- estimateMissSBM(A, Q, "node",       control = ctrl)$models[[1]]$ICL
  ICL_block_node <- estimateMissSBM(A, Q, "block-node", control = ctrl)$models[[1]]$ICL

  cl_A <- apply(stA$Z, 1, which.max)
  do.call(rbind, lapply(c("A", "B", "C"), function(m) {
    st <- get(paste0("st", m))
    clm <- apply(st$Z, 1, which.max)
    p <- environment(st$samp$update_parameters)$private$psi
    data.frame(g, n = n, method = m, observed_nodes = mean(pon$observedNodes),
               ARI = aricode::ARI(clm, s0$memberships), ARI_vs_A = aricode::ARI(clm, cl_A),
               rmse_pi = rmse_pairs(clm, st$theta, s0$memberships),
               a_hat = p[1], b_hat = p[2], a = psi[1], b = psi[2],
               J = J_zeta(st), ICL = ICL_degree(st),
               ICL_node = ICL_node, ICL_block_node = ICL_block_node,
               iterations = length(st$trace) - 1, decreases = st$decreases,
               time = unname(c(A = t_A, B = t_B, C = t_C)[paste0(m, ".elapsed")]))
  }))
}

t0  <- Sys.time()
res <- do.call(rbind, mclapply(seq_len(nrow(grid)), run, mc.cores = cores))
attr(res, "commit") <- system("git rev-parse --short HEAD", intern = TRUE)
dir.create("inst/simus/issue55-degree/results", showWarnings = FALSE)
saveRDS(res, sprintf("inst/simus/issue55-degree/results/results-n%d.rds", n))
cat(nrow(res) / 3, "networks,", format(round(difftime(Sys.time(), t0, units = "mins"), 1)), "\n")
