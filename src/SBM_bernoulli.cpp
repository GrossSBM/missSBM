#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

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

// M-step for the covariate Bernoulli model: maximizes, over (Gamma, beta), the weighted
// logistic log-likelihood:
//   J(Gamma, beta) = sum_{(i,j) observed} sum_{q,l} w_ijql * (Y_ij*eta_ijql - log(1+exp(eta_ijql)))
// with eta_ijql = Gamma_ql + beta'X_ij and w_ijql = Z_iq * Z_jl.
//
// The logit link being canonical for the Bernoulli family, J is globally concave in
// (Gamma, beta) jointly, so Newton-Raphson (eq. IRLS) converges reliably and quadratically.
//
// The Hessian has an "arrowhead" structure: d2J/dGamma_ql^2 is diagonal (no coupling between
// different (q,l) pairs), bordered by a Q^2 x K (or reduced x K for the symmetric case) coupling
// block and a dense K x K block for beta. The diagonal Gamma block lets each Newton step be
// solved via a Schur complement on a K x K system (K = number of covariates, typically small)
// instead of inverting the full (Q^2+K) x (Q^2+K) system.
//
// Symmetric (undirected) case: reduce to the actual Q(Q+1)/2 free parameters (q <= l),
// summing the (q,l) and (l,q) contributions onto the shared parameter, giving an exact
// Newton step on the true parameter space and exact symmetry preservation at every iteration.

// [[Rcpp::export]]
Rcpp::List M_step_sparse_bernoulli_covariates (
    Rcpp::List init_param,
    const arma::sp_mat& Y,
    const arma::sp_mat& R,
    const arma::cube&   X,
    const arma::mat&    Z,
    const bool sym,
    const int    maxIter = 50,
    const double tol     = 1e-10) {

    arma::mat Gamma = Rcpp::as<arma::mat>(init_param["Gamma"]); // (Q,Q)
    arma::vec beta  = Rcpp::as<arma::vec>(init_param["beta"]);  // (K,1)

    uword Q = Z.n_cols;
    uword K = X.n_slices;

    // index map for the (reduced, if sym) set of free Gamma parameters
    uword P = sym ? (Q * (Q + 1)) / 2 : Q * Q;
    arma::umat pair_ql(P, 2);
    {
      uword p = 0;
      for (uword q = 0; q < Q; q++) {
        for (uword l = (sym ? q : 0); l < Q; l++) {
          pair_ql(p, 0) = q; pair_ql(p, 1) = l; p++;
        }
      }
    }

    auto eval_loglik = [&](const arma::mat & Gm, const arma::vec & bt) -> double {
      double loglik = 0;
      sp_mat::const_iterator Rij = R.begin(), Rij_end = R.end();
      for (; Rij != Rij_end; ++Rij) {
        uword i = Rij.row(), j = Rij.col();
        arma::vec phi = X.tube(i, j);
        double mu  = as_scalar(bt.t() * phi);
        double yij = Y(i, j);
        for (uword q = 0; q < Q; q++) {
          for (uword l = 0; l < Q; l++) {
            double zz = Z(i, q) * Z(j, l);
            if (zz == 0.0) continue;
            double eta = Gm(q, l) + mu;
            loglik += zz * (yij * eta - std::log1p(std::exp(eta)));
          }
        }
      }
      return loglik;
    };

    // Levenberg-Marquardt-style damping: ridge grows across (rare) failed attempts within an
    // iteration -- e.g. when some block pair (q,l) has near-zero curvature, which routinely
    // happens once the VEM has nearly converged to hard clusters -- and cools back down once a
    // damped step succeeds, so the common case stays a cheap, undamped, quadratically
    // convergent Newton step.
    double ridge = 1e-10;
    double J_old = eval_loglik(Gamma, beta);
    int iterate = 0;
    bool converged = false;

    for (; iterate < maxIter; iterate++) {

      arma::vec g_shared(P, fill::zeros); // dJ/dGamma_shared(p)
      arma::vec d_shared(P, fill::zeros); // d2J/dGamma_shared(p)^2  (<= 0, undamped)
      arma::mat C(P, K, fill::zeros);     // d2J/dGamma_shared(p) dbeta
      arma::vec g_beta(K, fill::zeros);
      arma::mat H_beta(K, K, fill::zeros); // undamped

      sp_mat::const_iterator Rij = R.begin(), Rij_end = R.end();
      for (; Rij != Rij_end; ++Rij) {
        uword i = Rij.row(), j = Rij.col();
        arma::vec phi = X.tube(i, j);
        double mu  = as_scalar(beta.t() * phi);
        double yij = Y(i, j);

        for (uword p = 0; p < P; p++) {
          uword q = pair_ql(p, 0), l = pair_ql(p, 1);
          double zz = (sym && q != l) ? (Z(i,q)*Z(j,l) + Z(i,l)*Z(j,q)) : Z(i,q)*Z(j,l);
          if (zz == 0.0) continue;
          double eta  = Gamma(q, l) + mu;
          double prob = 1.0 / (1.0 + std::exp(-eta));
          double d = zz * (yij - prob);
          double v = zz * prob * (1.0 - prob);

          g_shared(p) += d;
          d_shared(p) -= v;
          for (uword k = 0; k < K; k++) C(p, k) -= v * phi(k);
          g_beta += d * phi;
          H_beta -= v * (phi * phi.t());
        }
      }

      bool accepted = false;
      arma::mat Gamma_new; arma::vec beta_new; double J_new = J_old;
      for (int damp_try = 0; damp_try < 30 && !accepted; damp_try++) {

        arma::vec d_damped = d_shared - ridge;
        arma::mat H_beta_damped = H_beta;
        H_beta_damped.diag() -= ridge;

        // Newton step solves H*delta = -g ; via Schur complement on the (diagonal) Gamma block:
        arma::vec Dinv = 1.0 / d_damped;
        arma::mat CD = C.each_col() % Dinv;
        arma::mat Sc = H_beta_damped - C.t() * CD;
        arma::vec rhs = C.t() * (Dinv % g_shared) - g_beta;

        arma::vec delta_beta;
        if (!arma::solve(delta_beta, Sc, rhs)) {
          ridge *= 10;
          continue;
        }
        arma::vec delta_shared = -Dinv % (g_shared + C * delta_beta);

        arma::mat delta_Gamma(Q, Q, fill::zeros);
        for (uword p = 0; p < P; p++) {
          uword q = pair_ql(p, 0), l = pair_ql(p, 1);
          delta_Gamma(q, l) = delta_shared(p);
          if (sym) delta_Gamma(l, q) = delta_shared(p);
        }

        // backtracking line search: guarantees monotonic ascent for a given (damped) direction
        double alpha = 1.0;
        for (int ls_iter = 0; ls_iter < 20; ls_iter++) {
          Gamma_new = Gamma + alpha * delta_Gamma;
          beta_new  = beta  + alpha * delta_beta;
          J_new = eval_loglik(Gamma_new, beta_new);
          if (J_new > J_old - 1e-12) { accepted = true; break; }
          alpha *= 0.5;
        }
        if (!accepted) ridge *= 10; else ridge = std::max(ridge / 10, 1e-10);
      }
      if (!accepted) {
        throw Rcpp::exception("M_step_sparse_bernoulli_covariates: Newton solver failed "
                               "to find an ascent direction (degenerate Hessian).");
      }

      Gamma = Gamma_new;
      beta  = beta_new;

      if (std::abs(J_new - J_old) < tol * (std::abs(J_new) + 1e-8)) {
        J_old = J_new;
        iterate++;
        converged = true;
        break;
      }
      J_old = J_new;
    }

    return Rcpp::List::create(
        Rcpp::Named("status", converged ? 1 : 0),
        Rcpp::Named("iterations", iterate),
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

  // NB: a version hoisting log(1+exp(Gamma+M(i,j))) into one Q x Q evaluation per dyad
  // (instead of one per (dyad, q) pair) was benchmarked and found slower in practice for
  // realistic (small) Q -- see the note in vLL_complete_sparse_bernoulli_covariates above.
  uword Q = Z.n_cols ;
  sp_mat::const_iterator Rij     = R.begin();
  sp_mat::const_iterator Rij_end = R.end();

  arma::mat log_tau = Y * Z * Gamma.t() + Y.t() * Z * Gamma ;
  log_tau.each_row() += log(pi) ;

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
