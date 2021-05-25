#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(nloptr)]]
// [[Rcpp::plugins(cpp11)]]

#include "nlopt_wrapper.h"
#include "packing.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------
//
// lower bound of the expectation of the complete log-likelihood

// [[Rcpp::export]]
double vLL_complete_sparse_bernoulli_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z,
    const arma::mat& theta,
    const arma::vec& pi
) {
  return(
    accu(Z.t() * Y * Z % log(theta/(1-theta))) +
      accu(Z.t() * R * Z % log(1-theta)) + accu(Z * log(pi))
  ) ;
}

// [[Rcpp::export]]
double vLL_complete_sparse_bernoulli_covariates(
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
        loglik += (*Yij) * Z(Yij.row(),q) * Z(Yij.col(),l) * (Gamma(q,l) + M(Yij.row(),Yij.col()));
      }
    }
  }

  for(; Rij != Rij_end; ++Rij) {
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
Rcpp::List M_step_sparse_bernoulli_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat& Z,
    const bool sym = true
) {

  arma::mat ZtYZ = Z.t() * Y * Z;
  arma::mat ZtRZ = Z.t() * R * Z;
  arma::mat res;

  if (sym) {
    res = (ZtYZ + ZtYZ.t()) / (ZtRZ + ZtRZ.t()) ;
  } else {
    res = ZtYZ / ZtRZ ;
  }
  res.replace(arma::datum::nan, 0) ;
  res.replace(0,     exp(arma::datum::log_min)) ;
  res.replace(1, 1 - exp(arma::datum::log_min)) ;

  arma::rowvec pi = mean(Z,0) ;
  pi.replace(0,     exp(arma::datum::log_min)) ;
  pi.replace(1, 1 - exp(arma::datum::log_min)) ;

  return Rcpp::List::create(
    Rcpp::Named("theta", Rcpp::List::create(Rcpp::Named("mean", wrap(res)))),
    Rcpp::Named("pi"   , as<NumericVector>(wrap(pi)))
  );
}

// [[Rcpp::export]]
Rcpp::List M_step_sparse_bernoulli_covariates (
    Rcpp::List init_param,
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::cube&   X,
    const arma::mat&    Z,
    const bool sym,
    Rcpp::List configuration) {

     // Conversion from R, prepare optimization
    const auto init_Gamma = Rcpp::as<arma::mat>(init_param["Gamma"]); // (Q,Q)
    const auto init_beta  = Rcpp::as<arma::vec>(init_param["beta"]);  // (M,1)

    const auto metadata = tuple_metadata(init_Gamma, init_beta);
    enum { GAMMA_ID, BETA_ID }; // Names for metadata indexes

    auto parameters = std::vector<double>(metadata.packed_size);
    metadata.map<GAMMA_ID>(parameters.data()) = init_Gamma;
    metadata.map<BETA_ID>(parameters.data())  = init_beta;

    auto optimizer = new_nlopt_optimizer(configuration, parameters.size());

    // Optimize
    auto objective_and_grad = [&metadata, &Y, &R, &X, &Z, &sym](const double * params, double * grad) -> double {

      const arma::mat gamma = metadata.map<GAMMA_ID>(params);
      const arma::vec beta = metadata.map<BETA_ID>(params);

      uword Q = Z.n_cols;
      uword K = X.n_slices;
      double loglik = 0;

      arma::mat gr_gamma = zeros<mat>(Q,Q);
      arma::vec gr_beta  = zeros<vec>(K);

      sp_mat::const_iterator Rij     = R.begin();
      sp_mat::const_iterator Rij_end = R.end();

        for(; Rij != Rij_end; ++Rij) {
          arma::vec phi = X.tube(Rij.row(),Rij.col());
          double mu = as_scalar(beta.t() * phi) ;
          arma::mat ZiqZjl = Z.row(Rij.row()).t() * Z.row(Rij.col()) ;
          loglik += accu(ZiqZjl % ( Y(Rij.row(), Rij.col()) * (gamma + mu) - log(1 + exp(gamma + mu ))) ) ;
          arma::mat delta = ZiqZjl % (Y(Rij.row(), Rij.col()) - 1 / (1 + exp(-(gamma + mu)))) ;
          if (sym) {
            delta = .5 * (delta + delta.t()) ;
          }
          gr_gamma += delta ;
          gr_beta  += accu(delta) * phi ;
        }

        metadata.map<GAMMA_ID>(grad) = -gr_gamma;
        metadata.map<BETA_ID>(grad)  = -gr_beta;

        return (-loglik);
    };
    OptimizerResult result = minimize_objective_on_parameters(optimizer.get(), objective_and_grad, parameters);

    arma::mat Gamma = metadata.copy<GAMMA_ID>(parameters.data());
    arma::vec beta  = metadata.copy<BETA_ID>(parameters.data());

    return Rcpp::List::create(
        Rcpp::Named("status", static_cast<int>(result.status)),
        Rcpp::Named("iterations", result.nb_iterations),
        Rcpp::Named("theta", Rcpp::List::create(Rcpp::Named("mean", 1/(1 + exp(-Gamma))))),
        Rcpp::Named("pi"   , as<NumericVector>(wrap(mean(Z,0)))),
        Rcpp::Named("beta" , as<NumericVector>(wrap(beta)))
    );
}

// -----------------------------------------------------------------
//
// EXPECTATION STEP

// [[Rcpp::export]]
Rcpp::NumericMatrix E_step_sparse_bernoulli_nocovariate(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat&  Z,
    const arma::mat&  theta,
    const arma::rowvec&  pi,
    const bool rescale = true) {

  // I use trans
  mat log_tau = Y * Z * trans(log(theta/(1 - theta))) + R * Z * trans(log(1 - theta)) +
    Y.t() * Z * log(theta/(1 - theta)) +  R.t() * Z * log(1 - theta) ;
  log_tau.each_row() += log(pi) ;

  if (rescale) {
    log_tau.each_row( [](rowvec& x){
      x = exp(x - max(x)) / sum(exp(x - max(x))) ;
    }) ;
  }

  return Rcpp::wrap(log_tau) ;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix E_step_sparse_bernoulli_covariates(
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::mat&  M,
    const arma::mat&  Z,
    const arma::mat&  Gamma,
    const arma::rowvec&  pi,
    const bool symmetric = true,
    const bool rescale   = true) {

  uword Q = Z.n_cols ;
  sp_mat::const_iterator Rij     = R.begin();
  sp_mat::const_iterator Rij_end = R.end();

  arma::mat log_tau = Y * Z * Gamma.t() + Y.t() * Z * Gamma ;
  log_tau.each_row() += log(pi) ;

  // constant in Tau(i_.) so useless
  //  log_tau.each_col() += sum( (Y % M) * Z, 1)  + sum( (Y % M).t() * Z, 1) ;

  for(; Rij != Rij_end; ++Rij) {
    for(arma::uword q=0; q < Q; q++){
       log_tau(Rij.row(), q) -= accu(Z.row(Rij.col()) % log (1 + exp(Gamma.row(q) + M(Rij.row(),Rij.col())))) ;
       if (symmetric) {
         log_tau(Rij.col(), q) -= accu(Z.row(Rij.row()) % log (1 + exp(Gamma.row(q) + M(Rij.col(),Rij.row())))) ;
       }
    }
  }

  if( rescale) {
    log_tau.each_row( [](arma::rowvec& x){
      x = exp(x - max(x)) / sum(exp(x - max(x))) ;
    }) ;
  }

  return Rcpp::wrap(log_tau) ;
}

