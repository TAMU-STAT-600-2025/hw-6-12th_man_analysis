LRMultiClass_HW3 <- function(X,
                             y,
                             beta_init = NULL,
                             numIter = 50,
                             eta = 0.1,
                             lambda = 1) {
  n <- nrow(X)
  p <- ncol(X)
  K <- length(unique(y))
  
  if (!all(X[, 1] == rep(1, n))) {
    stop("The first column of X should be 1s.")
  }
  if (n != length(y)) {
    stop("The number of rows in X should be equal to the length of y.")
  }
  if (eta <= 0) {
    stop("eta should be positive.")
  }
  if (lambda < 0) {
    stop("lambda should be non-ngeative.")
  }
  if (is.null(beta_init)) {
    beta_init <- matrix(0, nrow = p, ncol = K)
  } else {
    if (nrow(beta_init) != p) {
      stop("The number of rows in beta_init should be p.")
    }
    if (ncol(beta_init) != K) {
      stop("The number of columns in beta_init should be K.")
    }
  }
  
  beta <- beta_init
  objective <- vector(mode = "numeric", length = numIter + 1)
  
  Xbeta <- exp(X %*% beta)
  pk <- t(Xbeta / rowSums(Xbeta))
  
  idx <- matrix(0, n, K)
  idx[cbind(1:n, y + 1)] <- 1
  objective[1] <- -sum(idx * log(t(pk))) + (lambda / 2) * sum(beta ** 2)
  
  for (iter in 1:numIter) {
    for (k in 1:K) {
      wk <- pk[k, ] * (1 - pk[k, ])
      a <- crossprod(X * wk, X) + diag(lambda, p)
      b <- crossprod(X, pk[k, ] - ifelse(y == k - 1, 1, 0)) + lambda * beta[, k]
      beta[, k] <- beta[, k] - eta * solve(a, b)
    }
    
    Xbeta <- exp(X %*% beta)
    pk <- t(Xbeta / rowSums(Xbeta))
    
    objective[iter + 1] <- -sum(idx * log(t(pk))) + (lambda / 2) * sum(beta ** 2)
  }
  
  return(list(beta = beta, objective =  objective))
}

test_that("multiplication works", {
  # Generate mocking data
  X <- as.matrix(iris[1:4])
  y <- as.numeric(iris$Species) - 1
  n <- nrow(X)
  X <- cbind(rep(1, n), X)
  p <- ncol(X)
  K <- max(y) + 1
  beta_init <- matrix(0, nrow = p, ncol = K)
  
  # Test if the C++ implementation matches R's output
  expect_equal(LRMultiClass_HW3(X, y, beta_init),
               LRMultiClass(X, y, beta_init),
               ignore_attr = TRUE)
})
