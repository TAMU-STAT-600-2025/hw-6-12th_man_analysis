// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// For simplicity, no test data, only training data, and no error calculation.
// X - n x p data matrix
// y - n length vector of classes, from 0 to K-1
// numIter - number of iterations, default 50
// eta - damping parameter, default 0.1
// lambda - ridge parameter, default 1
// beta_init - p x K matrix of starting beta values (always supplied in right format)
// [[Rcpp::export]]
Rcpp::List LRMultiClass_c(const arma::mat& X, const arma::uvec& y,
                          const arma::mat& beta_init, int numIter = 50,
                          double eta = 0.1, double lambda = 1) {
  // All input is assumed to be correct
  
  // Initialize some parameters
  int K = max(y) + 1;  // number of classes
  int p = X.n_cols;
  int n = X.n_rows;
  arma::mat beta =
    beta_init;  // to store betas and be able to change them if needed
  arma::vec objective(numIter + 1);  // to store objective values
  
  // Initialize anything else that you may need
  arma::mat Xbeta = arma::exp(X * beta);
  arma::mat pk = (Xbeta.each_col() / arma::sum(Xbeta, 1)).t();
  arma::mat idx = arma::zeros<arma::mat>(n, K);
  for (int i = 0; i < n; i++) {
    idx(i, y(i)) = 1.;
  }
  objective[0] = -arma::accu(idx % arma::log(pk.t())) +
    (lambda / 2.) * arma::accu(arma::square(beta));
  
  // Newton's method cycle - implement the update EXACTLY numIter iterations
  for (int iter = 1; iter <= numIter; iter++) {
    // Update beta
    for (int k = 0; k < K; k++) {
      arma::vec pk_k = pk.row(k).t();
      arma::vec wk = pk_k % (1 - pk_k);
      arma::mat a = (X.each_col() % wk).t() * X + lambda * arma::eye(p, p);
      arma::mat b = X.t() * (pk_k - (y == k)) + lambda * beta.col(k);
      beta.col(k) -= eta * arma::solve(a, b);
    }
    // Update pk
    Xbeta = arma::exp(X * beta);
    pk = (Xbeta.each_col() / arma::sum(Xbeta, 1)).t();
    
    // Update objective value
    objective[iter] = -arma::accu(idx % arma::log(pk.t())) +
      (lambda / 2) * arma::accu(arma::square(beta));
  }
  
  // Create named list with betas and objective values
  return Rcpp::List::create(Rcpp::Named("beta") = beta,
                            Rcpp::Named("objective") = objective);
}
