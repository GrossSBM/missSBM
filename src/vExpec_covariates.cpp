#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double vExpec_covariates(IntegerMatrix Y, NumericMatrix gamma, arma::vec beta, arma::cube cov, NumericMatrix Tau, NumericVector alpha) {

  int N = Y.ncol();
  int Q = Tau.ncol();
  double loglik = 0;

  for(int i=0; i < N; i++) {
    for(int q=0; q < Q; q++){
      for (int j=0; j < N; j++) {
        for (int l=0; l < Q; l++) {
          if (j < i) {
            arma::vec param = cov.tube(i,j);
            loglik = loglik + arma::as_scalar(Tau(i,q)*Tau(j,l)*( (Y(i,j)-1)*(gamma(q,l)+beta.t()*param) + g(arma::as_scalar(gamma(q,l)+beta.t()*param)) ));
          }
        }
      }
      loglik = loglik + arma::as_scalar(Tau(i,q)*log(alpha[q]));
    }
  }
  return loglik ;
}
