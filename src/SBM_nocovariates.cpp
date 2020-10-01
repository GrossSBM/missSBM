#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::NumericMatrix E_step_nocovariate(
    const arma::mat& Y,
    const arma::mat& Y_bar,
    const arma::mat&  theta,
    const arma::mat&  Tau,
    const arma::vec&  alpha,
    const double& log_lambda,
    int fixPointIter) {

  mat Tau_new ;

  for (int iter=0; iter < fixPointIter; iter++) {
    Tau_new = Y * Tau * trans(log(theta)) + Y_bar * Tau * trans(log(1 - theta)) + log_lambda ;
  }

  return Rcpp::wrap(Tau_new);
}
