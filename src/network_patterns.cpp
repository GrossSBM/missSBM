// The patterns of the observed dyads, of the missing ones and of the edges of a dense network
// with NA, read in place, without any n x n temporary (#49): see .network_patterns() in R.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
using namespace Rcpp;

template <int RTYPE>
static bool is_missing(typename traits::storage_type<RTYPE>::type x) {
  return traits::is_na<RTYPE>(x);
}

template <int RTYPE>
static List patterns_impl(const Matrix<RTYPE>& A) {
  const int n = A.ncol();

  // directed if A is not symmetric (NA matching NA)
  bool directed = false;
  for (int j = 0; j < n && !directed; ++j) {
    for (int i = 0; i < j; ++i) {
      const bool na_ij = is_missing<RTYPE>(A(i, j)), na_ji = is_missing<RTYPE>(A(j, i));
      if (na_ij != na_ji || (!na_ij && A(i, j) != A(j, i))) { directed = true; break; }
    }
  }

  // the valid dyads: i != j if directed, i < j otherwise, column by column: the row indices
  // (0-based) and the column pointers of the compressed sparse column format (dgCMatrix)
  std::vector<int> obs_i, miss_i, edge_i;
  IntegerVector obs_p(n + 1), miss_p(n + 1), edge_p(n + 1);
  bool has_zero = false, has_one = false, binary = true;
  for (int j = 0; j < n; ++j) {
    const int imax = directed ? n : j;
    for (int i = 0; i < imax; ++i) {
      if (i == j) continue;
      if (is_missing<RTYPE>(A(i, j))) {
        miss_i.push_back(i);
      } else {
        obs_i.push_back(i);
        const double v = A(i, j);
        if (v == 1) {
          has_one = true; edge_i.push_back(i);
        } else if (v == 0) {
          has_zero = true;
        } else {
          binary = false;
        }
      }
    }
    obs_p[j + 1] = obs_i.size(); miss_p[j + 1] = miss_i.size(); edge_p[j + 1] = edge_i.size();
    // the diagonal is not a dyad, but its values are checked too (as by the former check)
    if (!is_missing<RTYPE>(A(j, j))) {
      const double v = A(j, j);
      if (v == 1) has_one = true; else if (v == 0) has_zero = true; else binary = false;
    }
  }
  return List::create(
    Named("R_i") = wrap(obs_i),  Named("R_p") = obs_p,
    Named("S_i") = wrap(miss_i), Named("S_p") = miss_p,
    Named("Y_i") = wrap(edge_i), Named("Y_p") = edge_p,
    Named("directed") = directed, Named("binary") = binary,
    Named("has_zero") = has_zero, Named("has_one") = has_one
  );
}

// [[Rcpp::export]]
List network_patterns_dense(SEXP A) {
  switch (TYPEOF(A)) {
    case REALSXP: return patterns_impl<REALSXP>(NumericMatrix(A));
    case INTSXP:  return patterns_impl<INTSXP>(IntegerMatrix(A));
    case LGLSXP:  return patterns_impl<LGLSXP>(LogicalMatrix(A));
    default: stop("the network must be a numeric, integer or logical matrix");
  }
}

// -----------------------------------------------------------------
// The update of the imputations of the degree sampling, in one pass over the missing dyads
// (#57): see degreeSampling_fit$update_imputation(). For a missing dyad (i, j),
//   nu_ij = logistic(logit(PI_ij) + eta_ij), eta_ij = C_ij + C_ji - b,
//   C_ij = hk_i (a b + b^2 (1 + D_i - nu_ij)),
// with nu_ji the imputation of (j, i) in a directed network (0 when observed), bounded to
// [eps, 1 - eps] as check_boundaries().

// [[Rcpp::export]]
NumericVector degree_imputation_cpp(const NumericVector& PI, const IntegerVector& i, const IntegerVector& j,
                                    const NumericVector& nu, const IntegerVector& transpose,
                                    const NumericVector& hk, const NumericVector& D,
                                    double a, double b, bool directed) {
  const double eps = std::numeric_limits<double>::epsilon();
  const double ab = a * b, b2 = b * b;
  const int m = PI.size();
  NumericVector out(m);
  for (int k = 0; k < m; ++k) {
    const int ik = i[k] - 1, jk = j[k] - 1;
    double nu_ji = nu[k];
    if (directed) nu_ji = (transpose[k] == NA_INTEGER) ? 0 : nu[transpose[k] - 1];
    const double eta = hk[ik] * (ab + b2 * (1 + D[ik] - nu[k])) + hk[jk] * (ab + b2 * (1 + D[jk] - nu_ji)) - b;
    const double p = PI[k];
    double x = 1 / (1 + std::exp(-(std::log(p / (1 - p)) + eta)));
    if (std::isnan(x)) x = eps;
    out[k] = std::min(std::max(x, eps), 1 - eps);
  }
  return out;
}
