// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec MyKmeans_c(const arma::mat& X, int K,
                            const arma::mat& M, int numIter = 100){
    // All input is assumed to be correct
    
    // Initialize some parameters
    int n = X.n_rows;
    int p = X.n_cols;
    arma::uvec Y(n); // to store cluster assignments
    
    // Initialize any additional parameters if needed
    arma::mat Mcopy = M;
    
    // For loop with kmeans algorithm
    for (int i = 0; i < numIter; ++i) {
      // Calculate Euclidean distances
      arma::vec Xnorm = arma::sum(arma::square(X), 1);
      
      arma::vec Mnorm = arma::sum(arma::square(Mcopy), 1);
      
      arma::mat XMproduct = X * Mcopy.t();
      
      arma::mat dists = arma::repmat(Xnorm, 1, K) + arma::repmat(Mnorm.t(), n, 1)
        - 2.0 * (XMproduct);
      
      // Assign clusters
      for (int i = 0; i < n; ++i) {
        Y(i) = dists.row(i).index_min();
      }
      
      // Check if any cluster disappeard
      arma::uvec cluster_counts(K, arma::fill::zeros);
      for (int i = 0; i < n; ++i){
        cluster_counts(Y(i))++;
      }
      if (arma::any(cluster_counts == 0u)) {
        Rcpp::stop("Error: One of the clusters has disappeared.");
      }
      
      // Recompute centers
      arma::mat newM(K, p, arma::fill::zeros);
      
      for (int i = 0; i < n; ++i) {
        newM.row(Y(i)) += X.row(i);
      }
      
      for (int k = 0; k < K; ++k) {
        newM.row(k) /= static_cast<double>(cluster_counts(k));
      }
      
      // Check convergence
      if (arma::approx_equal(newM, Mcopy, "absdiff", 1e-12)) {
        Rcpp::Rcout << "Converged after " << (i + 1) << "iterations.\n";
        for (int i = 0; i < n; ++i) Y(i) += 1;
        return(Y);
      }
      
      Mcopy = newM;
    }
    
    Rcpp::Rcout << "Reached maximum number of iterations (" << numIter << ").\n";
    for (int i = 0; i < n; ++i) Y(i) += 1;
    
    // Returns the vector of cluster assignments
    return(Y);
}

