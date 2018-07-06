#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double objective_Mstep_covariates(arma::vec param, IntegerMatrix Y, arma::cube cov, NumericMatrix Tau, bool directed) {

  int N = Y.ncol();
  int Q = Tau.ncol();
  int M = cov.n_slices;
  double loglik = 0;

  arma::mat gamma(&param[0], Q, Q, true);
  arma::vec beta(&param[std::pow(Q,2)], M, true);

  if (directed) {
    for(int i=0; i < N; i++) {
      for(int q=0; q < Q; q++){
        for (int j=0; j < N; j++) {
          for (int l=0; l < Q; l++) {
            if (i != j) {
              arma::vec param = cov.tube(i,j);
              double rp = arma::as_scalar(beta.t()*param);
              loglik += + arma::as_scalar(Tau(i,q)*Tau(j,l)*( (Y(i,j)-1)*(gamma(q,l) + rp) +
                g(arma::as_scalar(gamma(q,l) + rp)) ));
            }
          }
        }
      }
    }
  } else {
    for(int i=0; i < N; i++) {
      for(int q=0; q < Q; q++){
        for (int j=0; j < N; j++) {
          for (int l=0; l < Q; l++) {
            if (j < i) {
              arma::vec param = cov.tube(i,j);
              double rp = arma::as_scalar(beta.t()*param);
              loglik += arma::as_scalar(Tau(i,q)*Tau(j,l)*( (Y(i,j)-1)*(gamma(q,l) + rp) +
                g(arma::as_scalar(gamma(q,l) + rp)) ));
            }
          }
        }
      }
    }
  }
  return loglik;
}

// [[Rcpp::export]]
NumericVector gradient_Mstep_covariates(arma::vec param, IntegerMatrix Y, arma::cube cov, NumericMatrix Tau, bool directed) {

  int N = Y.ncol();
  int Q = Tau.ncol();
  int M = cov.n_slices;
  arma::mat gradGamma(Q,Q);
  arma::vec gradBeta(M);
  double acc;

  arma::mat gamma(&param[0], Q, Q, true);
  arma::vec beta(&param[std::pow(Q,2)], M, true);

  if(directed) {
    for(int q=0; q<Q; q++) {
      for(int l=0; l<Q; l++){
        acc = 0;
        for (int i=0; i<N; i++) {
          for (int j=0; j<N; j++) {
            if (j!=i) {
              arma::vec param = cov.tube(i,j);
              double rp = arma::as_scalar(beta.t()*param);
              acc = acc + arma::as_scalar(Tau(i,q)*Tau(j,l)*( Y(i,j) - 1 + u(rp + gamma(q,l)) ));
            }
          }
        }
        gradGamma(q,l) = acc;
      }
    }
  } else {
    for(int q=0; q<Q; q++) {
      for(int l=0; l<=q; l++){
        acc = 0;
        for (int i=0; i<N; i++) {
          for (int j=0; j<i; j++) {
            arma::vec param = cov.tube(i,j);
            double rp = arma::as_scalar(beta.t()*param);
            acc = acc + arma::as_scalar(Tau(i,q)*Tau(j,l)*( Y(i,j) - 1 + u(rp + gamma(q,l)) ));
          }
        }
        gradGamma(q,l) = acc;
        gradGamma(l,q) = acc;
      }
    }
  }

  if(directed) {
    for(int k=0; k<M; k++){
      acc = 0;
      for(int q=0; q<Q; q++) {
        for(int l=0; l<Q; l++){
          for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
              if (j != i) {
                arma::vec param = cov.tube(i,j);
                double rp = arma::as_scalar(beta.t()*param);
                acc = acc + arma::as_scalar(Tau(i,q)*Tau(j,l)*( (Y(i,j) - 1 + u(rp + gamma(q,l)))*param[k] ));
              }
            }
          }
        }
      }
      gradBeta[k] = acc;
    }
  } else {
    for(int k=0; k<M; k++){
      acc = 0;
      for(int q=0; q<Q; q++) {
        for(int l=0; l<Q; l++){
          for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
              if (j < i) {
                arma::vec param = cov.tube(i,j);
                double rp = arma::as_scalar(beta.t()*param);
                acc = acc + arma::as_scalar(Tau(i,q)*Tau(j,l)*( (Y(i,j) - 1 + u(rp + gamma(q,l)))*param[k] ));
              }
            }
          }
        }
      }
      gradBeta[k] = acc;
    }
  }

  arma::vec gradGammaVec = vectorise(gradGamma);
  arma::vec out = join_cols(gradGammaVec, gradBeta);
  return wrap(out);
}
