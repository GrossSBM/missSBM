// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// vLL_complete_sparse_bernoulli_nocovariate
double vLL_complete_sparse_bernoulli_nocovariate(const arma::sp_mat& Y, const arma::sp_mat& R, const arma::mat& Z, const arma::mat& theta, const arma::vec& pi);
RcppExport SEXP _missSBM_vLL_complete_sparse_bernoulli_nocovariate(SEXP YSEXP, SEXP RSEXP, SEXP ZSEXP, SEXP thetaSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(vLL_complete_sparse_bernoulli_nocovariate(Y, R, Z, theta, pi));
    return rcpp_result_gen;
END_RCPP
}
// vLL_complete_sparse_bernoulli_covariates
double vLL_complete_sparse_bernoulli_covariates(const arma::sp_mat& Y, const arma::sp_mat& R, const arma::mat& M, const arma::mat& Z, const arma::mat& Gamma, const arma::vec& pi);
RcppExport SEXP _missSBM_vLL_complete_sparse_bernoulli_covariates(SEXP YSEXP, SEXP RSEXP, SEXP MSEXP, SEXP ZSEXP, SEXP GammaSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(vLL_complete_sparse_bernoulli_covariates(Y, R, M, Z, Gamma, pi));
    return rcpp_result_gen;
END_RCPP
}
// M_step_sparse_bernoulli_nocovariate
Rcpp::List M_step_sparse_bernoulli_nocovariate(const arma::sp_mat& Y, const arma::sp_mat& R, const arma::mat& Z, const bool sym);
RcppExport SEXP _missSBM_M_step_sparse_bernoulli_nocovariate(SEXP YSEXP, SEXP RSEXP, SEXP ZSEXP, SEXP symSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const bool >::type sym(symSEXP);
    rcpp_result_gen = Rcpp::wrap(M_step_sparse_bernoulli_nocovariate(Y, R, Z, sym));
    return rcpp_result_gen;
END_RCPP
}
// M_step_sparse_bernoulli_covariates
Rcpp::List M_step_sparse_bernoulli_covariates(Rcpp::List init_param, const arma::sp_mat& Y, const arma::sp_mat& R, const arma::cube& X, const arma::mat& Z, Rcpp::List configuration);
RcppExport SEXP _missSBM_M_step_sparse_bernoulli_covariates(SEXP init_paramSEXP, SEXP YSEXP, SEXP RSEXP, SEXP XSEXP, SEXP ZSEXP, SEXP configurationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type init_param(init_paramSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type configuration(configurationSEXP);
    rcpp_result_gen = Rcpp::wrap(M_step_sparse_bernoulli_covariates(init_param, Y, R, X, Z, configuration));
    return rcpp_result_gen;
END_RCPP
}
// E_step_sparse_bernoulli_nocovariate
Rcpp::NumericMatrix E_step_sparse_bernoulli_nocovariate(const arma::sp_mat& Y, const arma::sp_mat& R, const arma::mat& Z, const arma::mat& theta, const arma::rowvec& pi, const bool rescale);
RcppExport SEXP _missSBM_E_step_sparse_bernoulli_nocovariate(SEXP YSEXP, SEXP RSEXP, SEXP ZSEXP, SEXP thetaSEXP, SEXP piSEXP, SEXP rescaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const bool >::type rescale(rescaleSEXP);
    rcpp_result_gen = Rcpp::wrap(E_step_sparse_bernoulli_nocovariate(Y, R, Z, theta, pi, rescale));
    return rcpp_result_gen;
END_RCPP
}
// E_step_sparse_bernoulli_covariates
Rcpp::NumericMatrix E_step_sparse_bernoulli_covariates(const arma::sp_mat& Y, const arma::sp_mat& R, const arma::mat& M, const arma::mat& Z, const arma::mat& Gamma, const arma::rowvec& pi, const bool symmetric, const bool rescale);
RcppExport SEXP _missSBM_E_step_sparse_bernoulli_covariates(SEXP YSEXP, SEXP RSEXP, SEXP MSEXP, SEXP ZSEXP, SEXP GammaSEXP, SEXP piSEXP, SEXP symmetricSEXP, SEXP rescaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::sp_mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type pi(piSEXP);
    Rcpp::traits::input_parameter< const bool >::type symmetric(symmetricSEXP);
    Rcpp::traits::input_parameter< const bool >::type rescale(rescaleSEXP);
    rcpp_result_gen = Rcpp::wrap(E_step_sparse_bernoulli_covariates(Y, R, M, Z, Gamma, pi, symmetric, rescale));
    return rcpp_result_gen;
END_RCPP
}
// cpp_test_nlopt
bool cpp_test_nlopt();
RcppExport SEXP _missSBM_cpp_test_nlopt() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cpp_test_nlopt());
    return rcpp_result_gen;
END_RCPP
}
// cpp_test_packing
bool cpp_test_packing();
RcppExport SEXP _missSBM_cpp_test_packing() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cpp_test_packing());
    return rcpp_result_gen;
END_RCPP
}
// roundProduct
Rcpp::NumericMatrix roundProduct(Rcpp::List covariates_list, arma::vec beta);
RcppExport SEXP _missSBM_roundProduct(SEXP covariates_listSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type covariates_list(covariates_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(roundProduct(covariates_list, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_missSBM_vLL_complete_sparse_bernoulli_nocovariate", (DL_FUNC) &_missSBM_vLL_complete_sparse_bernoulli_nocovariate, 5},
    {"_missSBM_vLL_complete_sparse_bernoulli_covariates", (DL_FUNC) &_missSBM_vLL_complete_sparse_bernoulli_covariates, 6},
    {"_missSBM_M_step_sparse_bernoulli_nocovariate", (DL_FUNC) &_missSBM_M_step_sparse_bernoulli_nocovariate, 4},
    {"_missSBM_M_step_sparse_bernoulli_covariates", (DL_FUNC) &_missSBM_M_step_sparse_bernoulli_covariates, 6},
    {"_missSBM_E_step_sparse_bernoulli_nocovariate", (DL_FUNC) &_missSBM_E_step_sparse_bernoulli_nocovariate, 6},
    {"_missSBM_E_step_sparse_bernoulli_covariates", (DL_FUNC) &_missSBM_E_step_sparse_bernoulli_covariates, 8},
    {"_missSBM_cpp_test_nlopt", (DL_FUNC) &_missSBM_cpp_test_nlopt, 0},
    {"_missSBM_cpp_test_packing", (DL_FUNC) &_missSBM_cpp_test_packing, 0},
    {"_missSBM_roundProduct", (DL_FUNC) &_missSBM_roundProduct, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_missSBM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
