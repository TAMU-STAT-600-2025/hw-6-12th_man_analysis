MyKmeansHW3 <- function(X, K, M = NULL, numIter = 100){
  #Ensure X is a matrix and count rows and columns of X
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  Y = c()
  
  # Check whether M is NULL or not. If NULL, initialize based on K randomly selected points from X. 
  # If not NULL, check for compatibility with X dimensions and K.
  if(is.null(M)){
    centroids = sample(n, K, replace = FALSE)
    M = X[centroids, , drop = FALSE]
  } else{
    if(!is.matrix(M)) stop("M must be a matrix")
    if(nrow(M) != K) stop("M must have K rows")
    if(ncol(M) != p) stop("M must have p columns")
    if(K > n) stop("K cannot be greater than n")
  }
  
  # Implement K-means algorithm. 
  # It should stop when either 
  # (i) the centroids don't change from one iteration to the next (exactly the same), or
  # (ii) the maximal number of iterations was reached, or
  # (iii) one of the clusters has disappeared after one of the iterations (in which case the error message is returned)
  for (i in 1:numIter) {
    # Assign each point to nearest centroid
    Xnorm = rowSums(X^2)                 
    Mnorm = rowSums(M^2)                 
    XMproduct = tcrossprod(X, M)
    Euclideandistances = outer(Xnorm, Mnorm, "+") - 2 * XMproduct
    
    # Assign clusters
    Y = max.col(-Euclideandistances)  
    # Check if any cluster disappeared
    cluster_counts = tabulate(Y, nbins = K)
    if (any(cluster_counts == 0)) {
      stop("Error: One of the clusters has disappeared.")
    }
    
    # Recompute centers
    newM = matrix(0, nrow = K, ncol = p)
    newM <- rowsum(X, Y) / as.vector(table(Y))
    
    # Check convergence
    if (all(newM == M)) {
      message("Converged after ", i, " iterations.")
      return(Y)
    }
    
    M = newM
  }
  
  message("Reached maximum number of iterations (", numIter, ").")
  
  # Return the vector of assignments Y
  return(Y)
}



test_that("kmeanscpp matches kmeansr", {
  X = matrix(c(0, 0.1, 5, 5.1), ncol = 1)
  M = matrix(c(1,2), ncol = 1)
  expect_equal(MyKmeansHW3(X, 2, M), MyKmeans_c(X, 2, M),
               ignore_attr = TRUE)
  
  X1 = matrix(rnorm(50, mean = 0, sd = 0.3), ncol = 2)
  X2 = matrix(rnorm(50, mean = 3, sd = 0.3), ncol = 2)
  X = rbind(X1, X2)
  M = matrix(c(1,2,3,4), ncol = 2)
  expect_equal(MyKmeansHW3(X, 2, M), MyKmeans_c(X, 2, M),
               ignore_attr = TRUE)
}
          )
