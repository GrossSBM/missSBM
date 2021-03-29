#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include "utils.h"

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------
//
// lower bound of the expectation of the complete log-likelihood

// [[Rcpp::export]]
double vLL_complete_sparse_bernoulli_undirected_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z,
    const arma::mat&  theta,
    const arma::vec&  pi
  ) {
  return(
      .5 * accu(Z.t() * Y * Z % log(theta/(1-theta))) +
       .5 * accu(Z.t() * R * Z % log(1-theta)) + accu(Z * log(pi))
  ) ;
}

// [[Rcpp::export]]
double vLL_complete_sparse_bernoulli_directed_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z,
    const arma::mat&  theta,
    const arma::vec&  pi
  ) {
  return(
      accu(Z.t() * Y * Z % log(theta/(1-theta))) +
        accu(Z.t() * R * Z % log(1-theta)) + accu(Z * log(pi))
  ) ;
}

// [[Rcpp::export]]
double vLL_complete_sparse_bernoulli_undirected_covariates(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& M,
    const arma::mat& Z,
    const arma::mat& Gamma,
    const arma::vec& pi
  ) {

  uword Q = Z.n_cols ;
  double loglik = accu(Z * log(pi)) ;

  sp_mat::const_iterator Yij     = Y.begin();
  sp_mat::const_iterator Yij_end = Y.end();
  sp_mat::const_iterator Rij     = R.begin();
  sp_mat::const_iterator Rij_end = R.end();

  for(; Yij != Yij_end; ++Yij) {
    if (Yij.row()>Yij.col()) {
      for(uword q=0; q < Q; q++){
        for (uword l=0; l < Q; l++) {
          loglik += Z(Yij.row(),q) * Z(Yij.col(),l) * (Gamma(q,l) + M(Yij.row(),Yij.col()));
        }
      }
    }
  }

  for(; Rij != Rij_end; ++Rij) {
    if (Rij.row()>Rij.col()) {
      for(uword q=0; q < Q; q++){
        for (uword l=0; l < Q; l++) {
          loglik += - Z(Rij.row(),q) * Z(Rij.col(),l) * log(1 + exp(Gamma(q,l) + M(Rij.row(),Rij.col())));
        }
      }
    }
  }

  return loglik ;
}

// [[Rcpp::export]]
double vLL_complete_sparse_bernoulli_directed_covariates(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& M,
    const arma::mat& Z,
    const arma::mat& Gamma,
    const arma::vec& pi
  ) {

  uword Q = Z.n_cols ;
  double loglik = accu(Z * log(pi)) ;

  sp_mat::const_iterator Yij     = Y.begin();
  sp_mat::const_iterator Yij_end = Y.end();
  sp_mat::const_iterator Rij     = R.begin();
  sp_mat::const_iterator Rij_end = R.end();

  for(; Yij != Yij_end; ++Yij) {
    for(uword q=0; q < Q; q++){
      for (uword l=0; l < Q; l++) {
        loglik += Z(Yij.row(),q) * Z(Yij.col(),l) * (Gamma(q,l) + M(Yij.row(),Yij.col()));
        }
      }
  }

  for(; Rij != Rij_end; ++Rij)
   {
    for(uword q=0; q < Q; q++){
      for (uword l=0; l < Q; l++) {
        loglik += - Z(Rij.row(),q) * Z(Rij.col(),l) * log(1 + exp(Gamma(q,l) + M(Rij.row(),Rij.col())));
        }
      }
  }

  return loglik ;
}

// -----------------------------------------------------------------
//
// MAXIMIZATION STEP

// [[Rcpp::export]]
Rcpp::NumericMatrix M_step_sparse_bernoulli_undirected_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z
    ) {
  return(Rcpp::wrap(
      (Z.t() * Y * Z) / (Z.t() * R * Z) )
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix M_step_sparse_bernoulli_directed_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z
    ) {
  return(Rcpp::wrap(
      (Z.t() * Y * Z) / (Z.t() * R * Z) )
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix M_step_sparse_bernoulli_undirected_covariates(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z
    ) {
  return(Rcpp::wrap(
      (Z.t() * Y * Z) / (Z.t() * R * Z) )
  );
}

// [[Rcpp::export]]
Rcpp::NumericMatrix M_step_sparse_bernoulli_directed_covariates(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z
    ) {
  return(Rcpp::wrap(
      (Z.t() * Y * Z) / (Z.t() * R * Z) )
  );
}

// -----------------------------------------------------------------
//
// EXPECTATION STEP

// [[Rcpp::export]]
Rcpp::NumericMatrix E_step_sparse_bernoulli_undirected_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat&  Z,
    const arma::mat&  theta,
    const arma::rowvec&  pi,
    const double& log_lambda = 0) {

    arma::mat log_tau = Y * Z * log(theta/(1 - theta)) + R * Z * log(1 - theta) + log_lambda ;
    log_tau.each_row() += log(pi) ;
    log_tau.each_row( [](arma::rowvec& x){
       x = exp(x - max(x)) / sum(exp(x - max(x))) ;
    }) ;

    return Rcpp::wrap(log_tau) ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix E_step_sparse_bernoulli_directed_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat&  Z,
    const arma::mat&  theta,
    const arma::rowvec&  pi,
    const double& log_lambda = 0) {

   // I use trans
    mat log_tau = Y * Z * trans(log(theta/(1 - theta))) + R * Z * trans(log(1 - theta)) +
      Y.t() * Z * log(theta/(1 - theta)) +  R.t() * Z * log(1 - theta) + log_lambda ;
    log_tau.each_row() += log(pi) ;
    log_tau.each_row( [](rowvec& x){
       x = exp(x - max(x)) / sum(exp(x - max(x))) ;
    }) ;

    return Rcpp::wrap(log_tau) ;
}

