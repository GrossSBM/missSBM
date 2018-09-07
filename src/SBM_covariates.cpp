#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

#include "utils.h"

using namespace Rcpp;
using namespace arma;

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix roundProduct(arma::cube X, arma::vec beta) {

  int N = X.n_rows;
  arma::mat M = arma::zeros<arma::mat>(N,N);

  for (int k = 0; k < beta.size(); k++) {
    M += X.slice(k) * beta[k];
  }

  return Rcpp::wrap(M);
}

// [[Rcpp::export]]
double vExpec_covariates(
    Rcpp::IntegerMatrix Y,
    Rcpp::NumericMatrix roundProd,
    Rcpp::NumericMatrix gamma,
    Rcpp::NumericMatrix Tau,
    Rcpp::NumericVector alpha
  ) {

  int N = Y.ncol();
  int Q = Tau.ncol();
  double loglik = 0;

  for(int i=0; i < N; i++) {
    for(int q=0; q < Q; q++){
      for (int j=0; j < N; j++) {
        for (int l=0; l < Q; l++) {
          if (j < i) {
            loglik += as_scalar(Tau(i,q)*Tau(j,l)*( (Y(i,j)-1)*(gamma(q,l) + roundProd(i,j)) + g(gamma(q,l) + roundProd(i,j)) ) );
          }
        }
      }
      loglik = loglik + Tau(i,q) * std::log(alpha[q]);
    }
  }
  return loglik ;
}

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

// [[Rcpp::export]]
List Mstep_covariates_undirected(arma::vec param, IntegerMatrix Y, arma::cube cov, NumericMatrix Tau) {

  int N = Y.ncol();
  int Q = Tau.ncol();
  int M = cov.n_slices;
  double loglik = 0;

  arma::mat gamma(&param[0], Q, Q, true);
  arma::vec beta (&param[Q*Q], M, true);

  arma::mat gr_gamma = zeros<mat>(Q,Q) ;
  arma::vec gr_beta  = zeros<vec>(M);

  for(int q=0; q < Q; q++) {
      for (int l=0; l < Q; l++) {
          for(int i=0; i < N; i++) {
              for (int j=0; j < i; j++) {
                  arma::vec phi   = cov.tube(i,j);
                  double gamma_phi_beta = gamma(q,l) + as_scalar(beta.t() * phi);
                  loglik  += Tau(i,q)*Tau(j,l)*( (Y(i,j) - 1) * gamma_phi_beta +  g(gamma_phi_beta) );
                  gr_beta += Tau(i,q)*Tau(j,l)*( (Y(i,j) - 1  + g_prime(gamma_phi_beta) ) * phi );
                  if (l <= q) {
                      gr_gamma(q,l) += Tau(i,q)*Tau(j,l)*(  Y(i,j) - 1  + g_prime(gamma_phi_beta) ) ;
                  }
              }
          }
          gr_gamma(l,q) = gr_gamma(q,l) ;
      }
  }

  arma::vec gr_gamma_v = vectorise(gr_gamma);
  arma::vec grad = join_cols(gr_gamma_v, gr_beta);

  return List::create(Named("objective") = - loglik, Named("gradient")  = - grad);
}

// [[Rcpp::export]]
List Mstep_covariates_directed(arma::vec param, IntegerMatrix Y, arma::cube cov, NumericMatrix Tau) {

  int N = Y.ncol();
  int Q = Tau.ncol();
  int M = cov.n_slices;
  double loglik = 0;

  arma::mat gamma(&param[0], Q, Q, true);
  arma::vec beta (&param[Q*Q], M, true);

  arma::mat gr_gamma = zeros<mat>(Q,Q) ;
  arma::vec gr_beta = zeros<vec>(M);

  for(int q=0; q < Q; q++) {
      for (int l=0; l < Q; l++) {
          for(int i=0; i < N; i++) {
              for (int j=0; j < N; j++) {
                 if (i != j) {
                    arma::vec phi   = cov.tube(i,j);
                    double gamma_phi_beta = gamma(q,l) + as_scalar(beta.t() * phi);
                    loglik        += Tau(i,q)*Tau(j,l)*( (Y(i,j) - 1) * gamma_phi_beta +  g(gamma_phi_beta) );
                    gr_gamma(q,l) += Tau(i,q)*Tau(j,l)*(  Y(i,j) - 1  + g_prime(gamma_phi_beta) ) ;
                    gr_beta       += Tau(i,q)*Tau(j,l)*( (Y(i,j) - 1  + g_prime(gamma_phi_beta) ) * phi );
                 }
              }
          }
          gr_gamma(l,q) = gr_gamma(q,l) ;
      }
  }

  arma::vec gr_gamma_v = vectorise(gr_gamma);
  arma::vec grad = join_cols(gr_gamma_v, gr_beta);

  return List::create(Named("objective") = - loglik, Named("gradient")  = - grad);
}
