# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  # Number of observations
  n <- nrow(Xtilde)
  # Compute residual sum of squares: ||Y - Xβ||²
  rss <- as.numeric(crossprod(Ytilde - Xtilde %*% beta))
  # Compute L1 penalty: sum of absolute coefficients
  l1_penalty <- sum(abs(beta))
  # Return LASSO objective: (1 / (2n)) * RSS + λ * ||β||₁
  return(rss / (2 * n) + lambda * l1_penalty)
}

#' LASSO with Proximal Gradient Descent
#'
#' @param X A n x p matrix of centered and scaled X
#' @param Y A n x 1 matrix of centered Y
#' @param lambda A tuning parameter
#' @param beta_start A p x 1 vector of starting point for coordinate-descent algorithm (default: NULL)
#' @param eps A precision level for convergence assessment (default: 0.0001)
#' @param s A step size for proximal gradient descent
#'
#' @returns A solution for LASSO regression and the optimal objective function value
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100 * 5), 100, 5)
#' beta <- c(-2, 0, 1, 0, 3)
#' epsilon <- rnorm(100)
#' Y <- X %*% beta + epsilon
#' beta_start <- rep(0, 5)
#' fitLASSO_prox_Nesterov(X, Y, 0.1, beta_start)
fitLASSO_prox_Nesterov <- function(X, Y, lambda, 
                                   beta_start = NULL, eps = 0.0001, s = 0.01){
  
  # Compatibility checks from ProximalExamples and initialization of beta_init
  p <- ncol(X)
  n <- nrow(X)
  if (length(Y) != n){
    stop("Dimensions of X and Y don't match")
  }
  if (lambda < 0){
    stop("Only non-negative values of lambda are allowed!")
  }
  if (is.null(beta_start)){
    # Initialize beta0
    beta_start <- rep(0,p)
  }else if (length(beta_start) != p){
    stop("Supplied initial starting point has length different from p!")
  }
  
  # Center and standardize X,Y as in HW4
  # Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  
  # Center X
  Xmeans <- colMeans(X)
  Xcentered <- X - matrix(Xmeans, n, p, byrow = TRUE)
  
  # Scale X
  weights = sqrt(colSums(Xcentered^2) / n)
  Xtilde <- Xcentered %*% diag(1 / weights)
  
  # Call C++ fitLASSOstandardized_prox_Nesterov_c function to implement the algorithm
  beta_tilde = fitLASSOstandardized_prox_Nesterov_c(Xtilde, Ytilde, lambda, beta_start, eps, s)
  
  fmin <- lasso(Xtilde, Ytilde, beta_tilde, lambda)
  
  # Perform back scaling and centering to get original intercept and coefficient vector
  
  beta <- beta_tilde / weights
  intercept <- Ymean - sum(Xmeans * beta)
  
  # Return 
  # beta - the solution (without center or scale)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin, intercept = intercept))
}