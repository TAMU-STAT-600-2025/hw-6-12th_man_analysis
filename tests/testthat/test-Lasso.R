library(Rcpp)
library(RcppArmadillo)

lasso_r <- function(Xtilde, Ytilde, beta, lambda){
  n <- nrow(Xtilde)
  rss <- as.numeric(crossprod(Ytilde - Xtilde %*% beta))
  l1_penalty <- sum(abs(beta))
  return(rss / (2 * n) + lambda * l1_penalty)
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