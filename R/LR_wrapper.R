#' Multi-Class Logistic Regression
#'
#' @param X n x p training data, 1st column should be 1s to account for intercept
#' @param y a vector of size n of class labels, from 0 to K-1
#' @param numIter a number of fixed iterations for the algorithm, default value is 50
#' @param eta learning rate, default value is 0.1
#' @param lambda ridge parameter, default value is 1
#' @param beta_init (optional) initial starting values of beta for the algorithm, should be p x K matrix
#'
#' @return p X K matrix of estimated beta values after numIter iterations
#' @export
#'
#' @examples
#' set.seed(123)
#' n <- 12
#' p <- 3
#' K <- 3
#' X_1 <- rnorm(n, mean = 0, sd = 1)
#' X_2 <- rnorm(n, mean = 0, sd = 0.5)
#' X <- cbind(1, X_1, X_2)
#' y <- sample(0:(K-1), n, replace = TRUE)
#' beta_hat <- LRMultiClass(X, y)
LRMultiClass <- function(X, y, beta_init = NULL, numIter = 50, eta = 0.1, lambda = 1){
  
  # Compatibility checks from HW3 and initialization of beta_init
  n <- nrow(X)
  p <- ncol(X)
  K <- max(y) + 1
  # Check that the first column of X and Xt are 1s, if not - display appropriate message and stop execution.
  # Check that first column of X are 1s
  if (!all(X[, 1] == 1)) {
    stop("The first column of X should be 1s.")
  }
  # Check for compatibility of dimensions between X and Y
  if (n != length(y)) {
    stop("The number of rows in X should be equal to the length of y.")
  }
  # Check eta is positive
  if (eta <= 0) {
    stop("eta should be positive.")
  }
  # Check lambda is non-negative
  if (lambda < 0) {
    stop("lambda should be non-negative")
  }
  
  # Check whether beta_init is NULL. If NULL, initialize beta with p x K matrix of zeroes. If not NULL, check for compatibility of dimensions with what has been already supplied.
  if (is.null(beta_init)) {
    beta_init = matrix(0, nrow = p, ncol = K)
  } else if (nrow(beta_init) != p) {
    stop("The number of rows in beta_init should be p.")
  } else if (ncol(beta_init) != K) {
    stop("The number of columns in beta_init should be K.")
  }
  
  # Call C++ LRMultiClass_c function to implement the algorithm
  out = LRMultiClass_c(X, y, beta_init, numIter, eta, lambda)
  
  # Return the class assignments
  return(out)
}