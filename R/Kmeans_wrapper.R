#' K-Means Clustering
#'
#' @param X Matrix of data points, where rows are observations (n) and columns are features (p)
#' @param K Integer specifying the number of clusters
#' @param M An optional K x p matrix of initial cluster centroids. If NULL (default), K points are randomly sampled from X as initial centroids
#' @param numIter The maximum number of iterations to run the algorithm (default: 100)
#'
#' @return Vector of length n containing the cluster assignments (integers from 1 to K) for each observation
#' @export
#'
#' @examples
#' X1 <- matrix(rnorm(50 * 2, mean = 0), 50, 2)
#' X2 <- matrix(rnorm(50 * 2, mean = 4), 50, 2)
#' data <- rbind(X1, X2)
#'
#' # Run K-Means with K=2
#' cluster_assignments <- MyKmeans(data, K = 2, numIter = 100)
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