#define ARMA_DONT_USE_OPENMP

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

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
arma::mat dist_l2(const arma::mat& M) {

  uword n = M.n_rows ;
  arma::mat D = trimatu(zeros(n, n)) ;
  for (uword i = 0; i < n; i++) {
    for (uword j = i + 1 ; j < n; j++) {
      D(i,j) = sqrt(accu(pow(M.row(i) - M.row(j), 2))) ;
    }
  }
  return(D) ;
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

// [[Rcpp::export]]
Rcpp::List spectral_clustering(const arma::sp_mat& A, const arma::vec& vBlocks) {

  uword n = A.n_cols ;

  // handling lonely souls
  arma::colvec d =  abs(A) * ones(A.n_rows, 1);
  // uvec unconnected = find(d == 0) ;
  // uvec connected   = find(d != 0) ;
  //
  // // remove unconnected nodes
  arma::sp_mat L = A;

  // compute normalized weighted Laplacian
  sp_mat::const_iterator Lij     = L.begin();
  sp_mat::const_iterator Lij_end = L.end();

  for(; Lij != Lij_end; ++Lij) {
    double normalize = d(Lij.row()) * d(Lij.col()) ;
    if (normalize > 0.0) {
      L(Lij.row(), Lij.col()) = -(*Lij) / sqrt(normalize);
    } else {
      L(Lij.row(), Lij.col()) = 0 ;
    }
  }

  // Normalized eigen values
  arma::vec eigval;
  arma::mat eigvec ;
  arma::eigs_sym(eigval, eigvec, L, vBlocks.max());

  Rcpp::List clustering(vBlocks.n_elem) ;

  for (uword k = 0; k < vBlocks.n_elem; k++) {
    if (vBlocks(k) == 1) {
      clustering[k] = arma::ones(n, 1) ;
    } else {

      arma::mat coordinates = eigvec.cols(0, vBlocks(k)-1) ;
      arma::mat dists = dist_l2(coordinates) ;
      uvec init_centroids = ind2sub(size(dists), dists.index_max()) ;

      while (init_centroids.n_elem < vBlocks(k)) {
        uword i = min(dists.rows(init_centroids)).index_max() ;
        init_centroids.resize(init_centroids.n_elem + 1) ;
        init_centroids.back() = i ;
      }

      arma::mat centroids = coordinates.rows(init_centroids) ;
      clustering[k] = kmeans_cpp(coordinates, centroids) + 1;
    }

  }

  return(clustering) ;
}

// [[Rcpp::export]]
arma::mat eigen_arma(const arma::sp_mat& L, const int& Kmax) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eigs_sym(eigval, eigvec, L, Kmax);
    return eigvec ;
}

