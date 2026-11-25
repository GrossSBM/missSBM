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

// Euclidean distance between rows i and j, accumulated as stats::dist() does
static double row_dist(const NumericMatrix& X, int i, int j) {
  double d = 0;
  for (int k = 0; k < X.ncol(); ++k) {
    const double dev = X(i, k) - X(j, k);
    d += dev * dev;
  }
  return std::sqrt(d);
}

// the farthest pair of rows, as the first maximum of as.matrix(dist(X)) in column-major order:
// (row, column), 1-based
// [[Rcpp::export]]
IntegerVector farthest_pair_cpp(const NumericMatrix& X) {
  const int n = X.nrow();
  double dmax = -1;
  int imax = 0, jmax = 0;
  // the first occurrence of the maximum in column-major order lies in the lower triangle
  for (int j = 0; j < n; ++j)
    for (int i = j + 1; i < n; ++i) {
      const double d = row_dist(X, i, j);
      if (d > dmax) { dmax = d; imax = i; jmax = j; }
    }
  return IntegerVector::create(imax + 1, jmax + 1);
}

// the distances of all the rows to row i (1-based), as a row of as.matrix(dist(X))
// [[Rcpp::export]]
NumericVector dist_to_row_cpp(const NumericMatrix& X, int i) {
  const int n = X.nrow();
  NumericVector d(n);
  for (int j = 0; j < n; ++j) d[j] = (j == i - 1) ? 0 : row_dist(X, j, i - 1);
  return d;
}
