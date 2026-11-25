#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
IntegerVector kmeans_cpp(const arma::mat & coordinates, arma::mat& input_centroids) {

    arma::mat centroids = input_centroids;

    arma::mat classif = zeros<mat>(coordinates.n_rows,1);

    arma::mat old_classif;

    unsigned int niter = 0;
    do
    {
        old_classif = classif;

        for(unsigned int i=0;i<classif.n_rows;i++)
        {

            double dmin = norm(coordinates.row(i)-centroids.row(0),2);
            unsigned int imin=0;
            for(unsigned int c=1;c<centroids.n_rows;c++)
            {
                double d = norm(coordinates.row(i)-centroids.row(c),2);
                if(d<dmin)
                {
                    dmin=d;
                    imin=c;
                }
            }


            classif(i) = imin;
        }

        arma::mat centroids_old = centroids;
        centroids.fill(0);

        colvec S(centroids.n_rows);
        S.fill(0);
        for(unsigned int i=0;i<classif.n_rows;i++)
        {
            centroids.row(classif(i)) += coordinates.row(i);
            S(classif(i))+=1;

        }

        for(unsigned int c=0;c<centroids.n_rows;c++)
        {
            // an empty cluster would otherwise divide by zero and propagate NaN;
            // keep its previous centroid instead
            if (S(c) > 0) {
                centroids.row(c) /= 1.0*S(c);
            } else {
                centroids.row(c) = centroids_old.row(c);
            }
        }

        niter++;

    } while(accu(classif-old_classif!=0)!=0 && niter<coordinates.n_rows);

    return Rcpp::wrap(classif);
}

// -----------------------------------------------------------------
// farthest-point seeding of kmeans_missSBM() without the n x n matrix of distances (#49)

// Euclidean distance between points i and j of Xt (d x n: one point per column, contiguous),
// accumulated over the coordinates in the same order as stats::dist()
static inline double point_dist(const double* Xt, int d, int i, int j) {
  const double* xi = Xt + (size_t) i * d;
  const double* xj = Xt + (size_t) j * d;
  double s = 0;
  for (int k = 0; k < d; ++k) {
    const double dev = xi[k] - xj[k];
    s += dev * dev;
  }
  return std::sqrt(s);
}

// the farthest pair of points, as the first maximum of as.matrix(dist(X)) in column-major order:
// (row, column), 1-based. Xt is t(X), d x n.
// [[Rcpp::export]]
IntegerVector farthest_pair_cpp(const NumericMatrix& Xt) {
  const int d = Xt.nrow(), n = Xt.ncol();
  const double* x = Xt.begin();
  double dmax = -1;
  int imax = 0, jmax = 0;
  // the first occurrence of the maximum in column-major order lies in the lower triangle
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) {
      const double dist = point_dist(x, d, i, j);
      if (dist > dmax) { dmax = dist; imax = i; jmax = j; }
    }
  return IntegerVector::create(imax + 1, jmax + 1);
}

// the distances of all the points to point i (1-based), as a row of as.matrix(dist(X)). Xt is
// t(X), d x n.
// [[Rcpp::export]]
NumericVector dist_to_row_cpp(const NumericMatrix& Xt, int i) {
  const int d = Xt.nrow(), n = Xt.ncol();
  const double* x = Xt.begin();
  NumericVector dist(n);
  for (int j = 0; j < n; ++j) dist[j] = (j == i - 1) ? 0 : point_dist(x, d, j, i - 1);
  return dist;
}
