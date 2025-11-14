#' Title
#'
#' @param X 
#' @param y 
#' @param numIter 
#' @param eta 
#' @param lambda 
#' @param beta_init 
#'
#' @return
#' @export
#'
#' @examples
#' # Give example
LRMultiClass <- function(X, y, beta_init = NULL, numIter = 50, eta = 0.1, lambda = 1){
  
  # Compatibility checks from HW3 and initialization of beta_init
  n <- nrow(X)
  p <- ncol(X)
  # Check that the first column of X and Xt are 1s, if not - display appropriate message and stop execution.
  # Check that first column of X are 1s
  if (!all(X[, 1] == 1)) {
    stop("X should have an intercept column of only ones")
  }
  # Check that first column of Xt are 1s
  if (!all(Xt[, 1] == 1)) {
    stop("Xt should have an intercept column of only ones")
  }
  # Check for compatibility of dimensions between X and Y
  if (n != length(y)) {
    stop("X and Y should have same number of rows")
  }
  # Check for compatibility of dimensions between Xt and Yt
  if (nrow(Xt) != length(yt)) {
    stop("Xt and Yt should have same number of rows")
  }
  # Check for compatibility of dimensions between X and Xt
  if (p != ncol(Xt)) {
    stop("X and Xt should have same number of columns")
  }
  # Check eta is positive
  if (eta <= 0) {
    stop("eta should be positive")
  }
  # Check lambda is non-negative
  if (lambda < 0) {
    stop("lambda should be non-negative")
  }
  
  # Check whether beta_init is NULL. If NULL, initialize beta with p x K matrix of zeroes. If not NULL, check for compatibility of dimensions with what has been already supplied.
  if (is.null(beta_init)) {
    beta_init = matrix(0, nrow = p, ncol = K)
  } else if (nrow(beta_init) != p) {
    stop("beta_init must have the same number of rows as columns of X")
  } else if (ncol(beta_init) != K) {
    stop("beta_init must have the same number of columns as the number of class labels in Y")
  }
  
  # Call C++ LRMultiClass_c function to implement the algorithm
  out = LRMultiClass_c(X, y, beta_init, numIter, eta, lambda)
  
  # Return the class assignments
  return(out)
}