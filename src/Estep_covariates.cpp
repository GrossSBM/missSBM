#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace Rcpp ;

// [[Rcpp::export]]
Rcpp::NumericMatrix E_step_covariates(
    Rcpp::IntegerMatrix Y,
    Rcpp::NumericMatrix roundProd,
    Rcpp::NumericMatrix gamma,
    Rcpp::NumericMatrix Tau,
    Rcpp::NumericVector alpha,
    int fixPointIter) {

  int N = Y.ncol();
  int Q = gamma.ncol();
  double acc;

  for (int iter=0; iter < fixPointIter; iter++) {
    for(int i=0; i < N; i++) {
      for(int q=0; q < Q; q++){

        acc = 0;
        for (int j=0; j < N; j++) {
          for (int l=0; l < Q; l++) {
            if (j != i) {
              acc = acc + Tau(j,l) * ( (Y(i,j) - 1) * ( gamma(q,l) + roundProd(i,j) ) + g( gamma(q,l) + roundProd(i,j) ) );
            }
          }
        }
        Tau(i,q) = alpha[q] * std::exp(acc);
      }
      Tau(i,_) = Tau(i,_)/sum(Tau(i,_));
    }
  }
  return Rcpp::wrap(Tau);
}

