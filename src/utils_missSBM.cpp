#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix roundProduct(Rcpp::List covariates_list, arma::vec beta) {

  uword N = Rcpp::as<mat>(covariates_list[0]).n_rows;
  uword P = Rcpp::as<mat>(covariates_list[0]).n_cols;
  arma::mat result = arma::zeros<arma::mat>(N,P);

  for (unsigned int k = 0; k < beta.size(); k++) {
    result += Rcpp::as<mat>(covariates_list[k]) * beta[k];
  }

  return Rcpp::wrap(result);
}

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

            double dmin = norm(rowvec(coordinates.row(i)-centroids.row(0)),2);
            unsigned int imin=0;
            for(unsigned int c=1;c<centroids.n_rows;c++)
            {
                double d = norm(rowvec(coordinates.row(i)-centroids.row(c)),2);
                if(d<dmin)
                {
                    dmin=d;
                    imin=c;
                }
            }


            classif(i) = imin;
        }

        centroids.fill(0);

        colvec S(classif.n_rows);
        S.fill(0);
        for(unsigned int i=0;i<classif.n_rows;i++)
        {
            centroids.row(classif(i)) += coordinates.row(i);
            S(classif(i))+=1;

        }

        for(unsigned int c=0;c<centroids.n_rows;c++)
        {
            centroids.row(c) /= 1.0*S(c);
        }

        niter++;

    } while(accu(classif-old_classif!=0)!=0 && niter<coordinates.n_rows);

    return Rcpp::wrap(classif);
}

Rcpp::IntegerVector k_means_with_arma(const arma::mat& coordinates, const uword & K) {

  uword n = coordinates.n_cols ;
  arma::mat centroids;
  bool status = arma::kmeans(centroids, coordinates, K, static_spread, 10, false);

  IntegerVector clustering(n) ;
  for(unsigned int i = 0; i < n; i++)
    {
      double dmin = norm(coordinates.col(i)-centroids.col(0),2);
      unsigned int imin=0;
      for(unsigned int c=1;c<centroids.n_cols;c++)
        {
          double d = norm(coordinates.col(i)-centroids.col(c),2);
          if(d<dmin)
            {
              dmin=d;
              imin=c;
            }
        }
      clustering(i) = imin + 1;
    }

  return(clustering);
}


//' Absolute Spectral Clustering
//'
//' @param A a dgCMatrix
//' @param vBlocks a vector of integer for the successive number of blocks considered
//'
//' @return a list of vector of clustering memberships
//'
//' @export
// [[Rcpp::export]]
Rcpp::List spectral_clustering_sparse(const arma::sp_mat& A, const arma::vec& vBlocks) {

  // wall_clock timer;
  // timer.tic();
  // double timing ;

  // initialization
  uword n = A.n_cols ;
  arma::sp_mat L = A ;
  if (!L.is_symmetric()) L = .5 * (L + L.t()) ;
  sp_mat::iterator Lij     = L.begin();
  sp_mat::iterator Lij_end = L.end();

  // timing = timer.toc();
  // std::cout << "init: number of seconds: " << timing << std::endl;

  // Compute normalized weighted Laplacian
  arma::colvec d =  abs(L) * ones(n, 1) ;
  for(; Lij != Lij_end; ++Lij) {
    *Lij = -(*Lij) / sqrt(d(Lij.row()) * d(Lij.col()));
  }

  // timing = timer.toc();
  // std::cout << "Laplacian: number of seconds: " << timing << std::endl;

  eigs_opts opts;
  opts.tol = 1e-4;

  // Normalized eigen values
  arma::vec eigval;
  arma::mat eigvec ;
  arma::eigs_sym(eigval, eigvec, L, vBlocks.max(), "lm", opts);

  // timing = timer.toc();
  // std::cout << "SVD: number of seconds: " << timing << std::endl;

  // k-means clustering for varying number of groups
  Rcpp::List clustering(vBlocks.n_elem) ;
  for (uword k = 0; k < vBlocks.n_elem; k++) {
    if (vBlocks(k) == 1) {
      IntegerVector ones(n) ;
      for(unsigned i = 0; i <n; ++i) { ones[i] = 1; }
      clustering[k] = ones ;
    } else {

      // we only consider K eigen vectors for K groups
      arma::mat coordinates = normalise(eigvec.cols(0, vBlocks(k)-1), 1, 1) ;

      clustering[k] = k_means_with_arma(coordinates.t(), vBlocks(k)) ;
    }

  }

  // timing = timer.toc();
  // std::cout << "kmeans: number of seconds: " << timing << std::endl;

  return(clustering) ;
}

//' Absolute Spectral Clustering
//'
//' @param A a matrix
//' @param vBlocks a vector of integer for the successive number of blocks considered
//'
//' @return a list of vector of clustering memberships
//'
//' @export
// [[Rcpp::export]]
Rcpp::List spectral_clustering_dense(const arma::mat& A, const arma::vec& vBlocks) {

  // initialization
  uword n = A.n_cols ;
  arma::mat L = A ;
  if (!L.is_symmetric()) L = .5 * (L + L.t()) ;

  // Compute normalized weighted Laplacian
  arma::colvec d =  1 / sqrt( abs(L) * ones(n, 1) ) ;
  L = - diagmat(d) * L * diagmat(d);

  // Normalized eigen values
  arma::vec eigval;
  arma::mat eigvec ;
  arma::eig_sym (eigval, eigvec, L);

  // k-means clustering for varying number of groups
  Rcpp::List clustering(vBlocks.n_elem) ;
  for (uword k = 0; k < vBlocks.n_elem; k++) {
    if (vBlocks(k) == 1) {
      IntegerVector ones(n) ;
      for(unsigned i = 0; i <n; ++i) { ones[i] = 1; }
      clustering[k] = ones ;
    } else {

      // we only consider K eigen vectors for K groups
      arma::mat coordinates = normalise(eigvec.cols(0, vBlocks(k)-1), 1, 1) ;

      clustering[k] = k_means_with_arma(coordinates.t(), vBlocks(k)) ;
    }

  }

  return(clustering) ;
}

