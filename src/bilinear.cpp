// Values of Z M Z' at the structural nonzeros of a sparse pattern (dgCMatrix), without the
// n x n product nor any (#nonzeros x Q) temporary.

#include <Rcpp.h>
using namespace Rcpp;

// ZMt = t(Z %*% M) and Zt = t(Z), Q x n (one node per column, contiguous); i and p: the row
// indices (0-based) and column pointers of the pattern. Returns, for each nonzero (i, j) in the
// order of the pattern, sum_q (Z M)_iq Z_jq.
// [[Rcpp::export]]
NumericVector bilinear_at_pattern_cpp(const NumericMatrix& ZMt, const NumericMatrix& Zt,
                                      const IntegerVector& i, const IntegerVector& p) {
  const int Q = Zt.nrow(), n = Zt.ncol();
  const double* zm = ZMt.begin();
  const double* z  = Zt.begin();
  NumericVector x(i.size());
  for (int j = 0; j < n; ++j) {
    const double* zj = z + (size_t) j * Q;
    for (int k = p[j]; k < p[j + 1]; ++k) {
      const double* zmi = zm + (size_t) i[k] * Q;
      double s = 0;
      for (int q = 0; q < Q; ++q) s += zmi[q] * zj[q];
      x[k] = s;
    }
  }
  return x;
}
