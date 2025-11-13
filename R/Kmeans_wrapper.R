#' Title
#'
#' @param X
#' @param K
#' @param M
#' @param numIter
#'
#' @return Explain return
#' @export
#'
#' @examples
#' # Give example
MyKmeans <- function(X, K, M = NULL, numIter = 100) {
  n = nrow(X) # number of rows in X
  
  # Check whether M is NULL or not.
  if (is.null(M)) {
    # If NULL, initialize based on K random points from X.
    M <- X[sample(1:nrow(X), K), ]
  } else {
    # If not NULL, check for compatibility with X dimensions.
    if (nrow(M) != K || ncol(M) != ncol(X)) {
      stop("Incompatible dimensions of M")
    }
    # Ensure M is a matrix
    M <- as.matrix(M)
  }
  
  # Call C++ MyKmeans_c function to implement the algorithm
  Y = MyKmeans_c(X, K, M, numIter)
  
  # Return the class assignments
  return(Y)
}