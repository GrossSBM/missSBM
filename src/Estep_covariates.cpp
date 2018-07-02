#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
NumericMatrix E_step_covariates(IntegerMatrix Y, arma::cube cov, NumericMatrix gamma, arma::vec beta, NumericMatrix Tau, NumericVector alpha) {

  int N = Y.ncol();
  int Q = gamma.ncol();
  double acc;

  for(int i=0; i < N; i++) {
    for(int q=0; q < Q; q++){
      acc = 0;
      for (int j=0; j < N; j++) {
        for (int l=0; l < Q; l++) {
          if (j != i) {
            arma::vec param = cov.tube(i,j);
            double rp = arma::as_scalar(beta.t()*param);
            acc = acc + arma::as_scalar(Tau(j,l)*((Y(i,j)-1)*(gamma(q,l)+rp) + g(arma::as_scalar(gamma(q,l)+rp))));
          }
        }
      }
      Tau(i,q) = alpha[q]*std::exp(acc);
    }
    Tau(i,_) = Tau(i,_)/sum(Tau(i,_));
  }
  return wrap(Tau);
}
