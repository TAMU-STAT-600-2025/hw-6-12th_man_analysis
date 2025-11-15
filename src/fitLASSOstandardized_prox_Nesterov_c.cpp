// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <cmath>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// v - p x 1 vector
// Soft-thresholding function, returns vector
// [[Rcpp::export]]
arma::colvec soft(const arma::colvec& v, double lambda){
  // Soft-threshold each coordinate (shrink by lambda, zero out small values)
  arma::colvec shrunken = arma::clamp(arma::abs(v) - lambda, 0, arma::datum::inf);
  return arma::sign(v) % shrunken;
}

// Xtilde - centered and scaled X, n x p
// Ytilde - centered Y, n x 1
// beta - value of beta at which to evaluate the function
// lamdba - tuning parameter
// Lasso objective function, returns scalar
// [[Rcpp::export]]
double lasso(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda){
  // Number of observations
  double n = Xtilde.n_rows;
  // Compute residuals: Y - Xβ
  arma::colvec residuals = Ytilde - Xtilde * beta;
  // Compute residual sum of squares: ||Y - Xβ||²
  double rss = arma::dot(residuals, residuals);
  // Compute L1 penalty: sum of absolute coefficients
  double l1 = arma::norm(beta, 1);
  // Return LASSO objective: (1 / (2n)) * RSS + λ * ||β||₁
  return(rss / (2.0 * n) + lambda * l1);
}

// lambda_t - scalar
// Update the Nesterov acceleration coefficient for the next iteration
// [[Rcpp::export]]
double calculate_lambda_tp1(double lambda_t){
  double discriminant = 1.0 + 4.0 * (lambda_t * lambda_t);
  return (1.0 + std::sqrt(discriminant)) / 2.0;
}


// Xtilde - centered and scaled X, n x p
// Ytilde - centered Y, n x 1
// lambda - tuning parameter
// beta0 - p vector of starting point for coordinate-descent algorithm, optional
// eps - precision level for convergence assessment, default 0.0001
// s - step size for proximal gradient
// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_prox_Nesterov_c(const arma::mat& Xtilde, const arma::colvec& Ytilde,
                                                  double lambda, const arma::colvec& beta_start, 
                                                double eps = 0.0001, double s = 0.01){
  // All input is assumed to be correct
  
  // Initialize some parameters
  int n = Xtilde.n_rows, p = Xtilde.n_cols;
  arma::colvec beta(p);
  
  beta = beta_start;
  double lambda_t = 1;
  arma::colvec xt_old(p);
  xt_old = beta;
  double f_old = lasso(Xtilde, Ytilde, beta, lambda);
  while (true) {
    arma::colvec gradient_step(p);
    gradient_step = beta + (s / n) * Xtilde.t() * (Ytilde - Xtilde * beta);
    arma::colvec xt(p);
    xt = soft(gradient_step, lambda * s);
    double lambda_tp1 = calculate_lambda_tp1(lambda_t);
    beta = xt + (lambda_t - 1) / (lambda_tp1) * (xt - xt_old);
    
    double f_new = lasso(Xtilde, Ytilde, beta, lambda);
    double error = f_old - f_new;
    if (std::abs(error) < eps) {
      break;
    }
    f_old = f_new;
    xt_old = xt;
    lambda_t = lambda_tp1;
  }

  return beta;
}
