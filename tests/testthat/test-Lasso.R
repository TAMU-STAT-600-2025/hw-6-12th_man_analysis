library(Rcpp)
library(RcppArmadillo)

lasso_r <- function(Xtilde, Ytilde, beta, lambda){
  n <- nrow(Xtilde)
  rss <- as.numeric(crossprod(Ytilde - Xtilde %*% beta))
  l1_penalty <- sum(abs(beta))
  return(rss / (2 * n) + lambda * l1_penalty)
}
standardizeXY <- function(X, Y){
  # [ToDo] Center Y
  # Compute mean of Y
  Ymean <- mean(Y)
  # Center Y by subtracting its mean
  Ytilde <- Y - Ymean
  # [ToDo] Center and scale X
  # Compute column means of X
  Xmeans <- colMeans(X)
  # Center X by subtracting its column means
  Xcentered <- sweep(X, 2, Xmeans, "-")
  # Compute weights that will be used in scaling X
  weights <- sqrt(colMeans(Xcentered * Xcentered))
  # Scale X so each column has unit variance
  Xtilde <- sweep(Xcentered, 2, weights, "/")
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  if (nrow(Xtilde) != length(Ytilde)){
    stop("Xtilde and Ytilde should have same number of rows")
  }
  #[ToDo]  Check that lambda is non-negative
  if (lambda < 0) {
    stop("lambda should be non-negative")
  }
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  
  if (is.null(beta_start)) {
    beta_start = rep(0, ncol(Xtilde))
  } else if (length(beta_start) != ncol(Xtilde)) {
    stop("beta_start must have the same number of entries as columns of Xtilde")
  }
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  
  # Initialize variables to track changes in the objective function across iterations
  f_new <- lasso(Xtilde, Ytilde, beta_start, lambda)
  # Ensure that first iteration of while loop runs
  f_old <- f_new + 2 * eps
  # Initialize beta
  beta <- beta_start
  # store number of observations
  n <- nrow(Xtilde)
  # store number of explanatory variables
  p <- ncol(Xtilde)
  # Check that difference between objective functions is greater than eps
  while (f_old - f_new >= eps){
    # Update previous objective value
    f_old <- f_new
    # Update previous value of beta
    beta_old <- beta
    XB <- Xtilde %*% beta_old
    # Calculate residual vector
    r <- Ytilde - XB
    for (j in 1:p){
      # gradient component for coordinate j: (1/n) * X_jᵀR 
      XjR_scaled <- sum(Xtilde[, j] * r) / n
      # Update β_j via soft-thresholding: S(β_old[j] + (X_jᵀ r)/n, λ)
      beta[j] <- soft(beta_old[j] + XjR_scaled, lambda)
      # Compute coefficient change for feature j: β_old[j] − β_new[j]
      delta_beta_j <- beta_old[j] - beta[j]
      # Update residuals: r ← r + X_j * (β_old[j] − β[j])
      r <- r + Xtilde[, j] * delta_beta_j
    }
    # Evaluate LASSO objective at the updated coefficients β
    f_new <- lasso(Xtilde, Ytilde, beta, lambda)
  }
  fmin = f_new
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  return(list(beta = beta, fmin = fmin))
}

fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  n <- nrow(Xtilde)
  if (n != length(Ytilde)){
    stop("Xtilde and Ytilde should have same number of rows")
  }
  p <- ncol(Xtilde)
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  
  if (!is.null(lambda_seq)){
    # Only keep lambdas >= 0
    lambda_seq <- lambda_seq[lambda_seq >= 0]
    # Sort lambdas from largest to smallest
    lambda_seq <- sort(lambda_seq, decreasing = TRUE)
    if(length(lambda_seq) == 0){
      warning("Invalid lambda_seq detected: only positive values should be supplied.")
      lambda_seq <- NULL
    }
  }
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution)
  if (is.null(lambda_seq)){
    XtY <- t(Xtilde) %*% Ytilde
    # Compute λ_max = max(|X_jᵀY| / n);
    # when λ ≥ λ_max, the soft-thresholding step sets all β_j = 0 (null solution)
    lambda_max <- max(abs(XtY)) / n
    # create a sequence of length n_lambda as
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  # calculate length of lambda_seq
  num_lambda <- length(lambda_seq)
  # initialize matrix of solutions at each lambda value
  beta_mat <- matrix(0, nrow = p, ncol = num_lambda)
  # initialize vector of objective function values at solution for each lambda
  fmin_vec <- rep(0, num_lambda)
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  beta <- rep(0, p)
  for (i in 1:num_lambda){
    lambda <- lambda_seq[i]
    # Perform coordinate descent LASSO for the current value of λ
    out <- fitLASSOstandardized(Xtilde, Ytilde, lambda = lambda, beta_start = beta, eps = eps)
    # Store the optimal coefficients and objective value for this λ
    beta_mat[, i] <- out$beta
    fmin_vec[i] <- out$fmin
    # Use the current solution as the warm start for the next λ
    beta <- out$beta
  }
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  scaledXY <- standardizeXY(X, Y)
  # Store standardized X, standardized Y, and the scaling weights used to give X unit variance
  Xtilde <- scaledXY$Xtilde
  Ytilde <- scaledXY$Ytilde
  weights <- scaledXY$weights
  # Store column means of untransformed X and Y
  Ymean <- scaledXY$Ymean
  Xmeans <- scaledXY$Xmeans
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  out <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, n_lambda = n_lambda, eps = eps)
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  # Rescale coefficient matrix back to the original data scale
  beta_mat <- sweep(out$beta_mat, 1, weights, "/")
  # β₀ = Ymean − Xmeansᵀ β
  beta0_vec <- Ymean - as.numeric(crossprod(Xmeans, beta_mat))
  lambda_seq <- out$lambda_seq
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


test_that("soft() correctly applies soft-thresholding to vectors", {
  a <- c(-3, -1, -0.5, 0.5, 1, 3)
  lambda <- 1
  out <- soft(a, lambda)
  
  expected <- c(-2, 0, 0, 0, 0, 2)
  
  expect_equal(out, expected, ignore_attr = TRUE)
})

test_that("lasso matches lasso_r for multiple inputs", {
  lambda <- 2
  
  # Example #1:
  set.seed(123)
  X1 <- matrix(rnorm(500), 100, 5)
  Y1 <- rep(0, 100)
  beta1 <- sample(1:10, 5)
  expect_equal(lasso_r(X1, Y1, beta1, lambda),
               lasso(X1, Y1, beta1, lambda),
               tolerance = 1e-8)
  
  # Example #2:
  set.seed(456)
  X2 <- matrix(rnorm(100), 50, 2)
  Y2 <- rep(50, 50)
  beta2 <- sample(-3:3, 2)
  expect_equal(lasso_r(X2, Y2, beta2, lambda),
               lasso(X2, Y2, beta2, lambda),
               tolerance = 1e-8)
  
  # Example #3:
  set.seed(789)
  X3 <- matrix(rnorm(1000), 10, 100)
  Y3 <- rnorm(10)
  beta3 <- -49:50  
  expect_equal(lasso_r(X3, Y3, beta3, lambda),
               lasso(X3, Y3, beta3, lambda),
               tolerance = 1e-8)
})

test_that("Check that calculate_lambda_tp1 returns correct output", {
  
  # --- Test Case 1: lambda_t = 0.0 ---
  # Expected: 1.0
  result <- calculate_lambda_tp1(0.0)
  expect_equal(result, 1.0, tolerance = 1e-6)
  
  
  # --- Test Case 2: lambda_t = sqrt(2) ---
  # Expected: 2.0
  result <- calculate_lambda_tp1(sqrt(2))
  expect_equal(result, 2.0, tolerance = 1e-6)
  
  
  # --- Test Case 3: lambda_t = sqrt(6) ---
  # Expected: 3.0
  result <- calculate_lambda_tp1(sqrt(6))
  expect_equal(result, 3.0, tolerance = 1e-6)
  
})

test_that("fitLASSOstandardized_prox_Nesterov_c matches R implementation", {
  
  lambda <- 0.1
  eps    <- 1e-14
  s      <- 0.01
  
  # Example 1:
  set.seed(101)
  Xtilde1 <- matrix(rnorm(40), nrow = 8, ncol = 5)
  Ytilde1 <- Xtilde1[, 2] + rnorm(8, sd = 10)
  centeredXY <- standardizeXY(Xtilde1, Ytilde1)
  cXtilde1 <- centeredXY$Xtilde
  cYtilde1 <- centeredXY$Ytilde
  beta_start1 <- rep(1, 5)
  
  expect_equal(
    as.vector(
      fitLASSOstandardized_prox_Nesterov_c(
        cXtilde1, cYtilde1, lambda, beta_start1, eps, s
      )
    ),
    fitLASSOstandardized(
      cXtilde1, cYtilde1, lambda, beta_start1, eps
    )$beta,
    tolerance = 1e-5
  )
  
  # Example 2:
  set.seed(202)
  Xtilde2 <- matrix(rnorm(72), nrow = 12, ncol = 6)
  Ytilde2 <- 2 * Xtilde2[, 1] - 0.5 * Xtilde2[, 4] + rnorm(12, sd = 5)
  centeredXY <- standardizeXY(Xtilde2, Ytilde2)
  cXtilde2 <- centeredXY$Xtilde
  cYtilde2 <- centeredXY$Ytilde
  beta_start2 <- rep(0, 6)
  
  expect_equal(
    as.vector(
      fitLASSOstandardized_prox_Nesterov_c(
        cXtilde2, cYtilde2, lambda, beta_start2, eps, s
      )
    ),
    fitLASSOstandardized(
      cXtilde2, cYtilde2, lambda, beta_start2, eps
    )$beta,
    tolerance = 1e-5
  )
  
  # Example 3:
  set.seed(303)
  Xtilde3 <- matrix(rnorm(90), nrow = 15, ncol = 6)
  Ytilde3 <- rowSums(Xtilde3) + rnorm(15, sd = 1)
  centeredXY <- standardizeXY(Xtilde3, Ytilde3)
  cXtilde3 <- centeredXY$Xtilde
  cYtilde3 <- centeredXY$Ytilde
  beta_start3 <- rep(-1, 6)
  
  expect_equal(
    as.vector(
      fitLASSOstandardized_prox_Nesterov_c(
        cXtilde3, cYtilde3, lambda, beta_start3, eps, s
      )
    ),
    fitLASSOstandardized(
      cXtilde3, cYtilde3, lambda, beta_start3, eps
    )$beta,
    tolerance = 1e-5
  )
  
})

test_that("fitLASSO_prox_Nesterov matches fitLASSO for a single lambda", {
  
  set.seed(123)
  n <- 50
  p <- 8
  
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  beta_true <- c(1.5, -2, 0.5, rep(0, p - 3))
  Y <- as.vector(X %*% beta_true + rnorm(n, sd = 1))
  
  lambda <- 0.1
  eps    <- 1e-12
  s      <- 0.01
  

  fit_path <- fitLASSO(
    X, Y,
    lambda_seq = lambda,
    n_lambda   = 1,
    eps        = eps
  )
  print(fit_path$beta_mat)
  
  beta_seq      <- fit_path$beta_mat
  intercept_seq <- as.numeric(fit_path$beta0_vec[1])
  

  fit_nest <- fitLASSO_prox_Nesterov(
    X, Y,
    lambda     = lambda,
    beta_start = rep(0, p),
    eps        = eps,
    s          = s
  )
  
  expect_equal(fit_nest$beta,      beta_seq,      tolerance = 1e-4)
  expect_equal(fit_nest$intercept, intercept_seq, tolerance = 1e-4)
})
