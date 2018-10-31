#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::NumericMatrix E_step_nocovariate(
    const arma::mat& Y,
    const arma::mat& Y_bar,
    const arma::mat&  pi,
    const arma::mat&  Tau,
    const arma::vec&  alpha,
    const double& log_lambda,
    int fixPointIter) {

  int N = Y.n_cols;
  int Q = alpha.n_elem;
  mat Tau_new ;
  double acc;

  for (int iter=0; iter < fixPointIter; iter++) {
    Tau_new = Y * Tau * trans(log(pi)) + Y_bar * Tau * trans(log(1 - pi)) + log_lambda ;
  }

  return Rcpp::wrap(Tau_new);
}
