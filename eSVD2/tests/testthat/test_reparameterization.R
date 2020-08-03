context("Testing reparameterization")

test_that(".identification works", {
 res <- .identification(diag(5), 2*diag(5))
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == 5))
})

test_that(".identification works on an instance where you need to worry about underflow", {
 load("../assets/identification1.RData")
 
 res <- .identification(cov_x, cov_y)
 expect_true(sum(abs(Im(res))) == 0)
})

test_that(".identification is correct for a deterministic setting", {
 cov_x <- matrix(c(2,1,1,2),2,2)
 cov_y <- matrix(c(5,-1,-1,5),2,2)
 res <- .identification(cov_x, cov_y)
 
 mat1 <- res %*% cov_x %*% t(res)
 mat2 <- solve(t(res)) %*% cov_y %*% solve(res)
 
 res <- .identification(diag(5), 2*diag(5))
 
 expect_true(sum(abs(mat1 - mat2)) <= 1e-6)
})

test_that(".identification issues warning for rank defficient setting", {
 set.seed(10)
 x_mat <- matrix(rnorm(12), ncol = 3, nrow = 4)
 y_mat <- matrix(rnorm(15), ncol = 3, nrow = 5)
 
 cov_x <- stats::cov(x_mat)
 cov_y <- stats::cov(y_mat)
 
 cov_y[,3] <- 0; cov_y[3,] <- 0
 expect_warning(.identification(cov_x, cov_y))
})

test_that(".identification is correct", {
 set.seed(20)
 cov_x <- cov(MASS::mvrnorm(n = 10, rep(0, 5), diag(5)))
 cov_y <- cov(MASS::mvrnorm(n = 10, rep(0, 5), toeplitz(5:1)))
 
 res <- .identification(cov_x, cov_y)
 
 cov_x2 <- res %*% cov_x %*% t(res)
 cov_y2 <- solve(t(res)) %*% cov_y %*% solve(res)
 
 expect_true(sum(abs(cov_x2 - cov_y2)) <= 1e-6)
 expect_true(abs(sum(cov_x2) - sum(diag(cov_x2))) <= 1e-6)
 expect_true(abs(sum(cov_y2) - sum(diag(cov_y2))) <= 1e-6)
})

test_that(".identification for 1-dim covariances (just variances)", {
 cov_x <- matrix(10, 1, 1)
 cov_y <- matrix(5, 1, 1)
 
 res <- .identification(cov_x, cov_y)
 
 cov_x2 <- res %*% cov_x %*% t(res)
 cov_y2 <- solve(t(res)) %*% cov_y %*% solve(res)
 
 expect_true(sum(abs(cov_x2 - cov_y2)) <= 1e-6)
})

####################

## .reparameterize is correct

test_that(".reparameterize works", {
 set.seed(10)
 x_mat <- MASS::mvrnorm(60, rep(0, 5), diag(5))
 y_mat <- MASS::mvrnorm(50, rep(1, 5), 2*diag(5))
 
 res <- .reparameterize(x_mat, y_mat)
 
 expect_true(is.list(res))
 expect_true(length(res) == 2)
 expect_true(all(dim(res$x_mat) == dim(x_mat)))
 expect_true(all(dim(res$y_mat) == dim(y_mat)))
})

test_that(".reparameterize works for rank 1 matrices", {
 set.seed(10)
 x_mat <- matrix(rnorm(50), ncol = 1)
 y_mat <-matrix(rnorm(50), ncol = 1)
 
 res <- .reparameterize(x_mat, y_mat)
 
 expect_true(length(res) == 2)
})

test_that(".reparameterize preserves the inner products", {
 trials <- 100
 
 bool_vec <- sapply(1:trials, function(x){
  set.seed(10*x)
  x_mat <- MASS::mvrnorm(10, rep(0, 5), diag(5))
  y_mat <- MASS::mvrnorm(10, rep(1, 5), 2*diag(5))
  
  res <- .reparameterize(x_mat, y_mat)
  
  pred_mat <- x_mat %*% t(y_mat)
  pred_mat2 <- res$x_mat %*% t(res$y_mat)
  
  sum(abs(pred_mat - pred_mat2)) <= 1e-6
 })
 
 expect_true(all(bool_vec))
})

test_that(".reparameterize yields the same second moment matrix", {
 trials <- 100
 
 bool_vec <- sapply(1:trials, function(x){
  set.seed(11*x)
  x_mat <- MASS::mvrnorm(5, rep(0, 5), diag(5))
  y_mat <- MASS::mvrnorm(5, rep(1, 5), 2*diag(5))
  
  res <- .reparameterize(x_mat, y_mat)
  
  cov_u <- t(res$x_mat) %*% res$x_mat
  cov_v <- t(res$y_mat) %*% res$y_mat
  
  sum(abs(cov_u - cov_v)) < 1e-6
 })
 
 expect_true(all(bool_vec))
})


test_that(".reparameterize yields the same covariance matrix", {
 trials <- 100
 
 bool_vec <- sapply(1:trials, function(x){
  set.seed(12*x)
  n <- 10
  p <- 20
  x_mat <- MASS::mvrnorm(n, rep(0, 5), diag(5))
  y_mat <- MASS::mvrnorm(p, rep(1, 5), 2*diag(5))
  
  res <- .reparameterize(x_mat, y_mat, equal_covariance = T)
  
  cov_u <- t(res$x_mat) %*% res$x_mat/n
  cov_v <- t(res$y_mat) %*% res$y_mat/p
  
  sum(abs(cov_u - cov_v)) < 1e-6
 })
 
 expect_true(all(bool_vec))
})

test_that(".reparameterize yields diagonal covariances", {
 trials <- 100
 
 bool_vec <- sapply(1:trials, function(x){
  set.seed(13*x)
  x_mat <- MASS::mvrnorm(10, rep(0, 5), diag(5))
  y_mat <- MASS::mvrnorm(20, rep(1, 5), 2*diag(5))
  
  res <- .reparameterize(x_mat, y_mat)
  
  cov_u <- t(res$x_mat) %*% res$x_mat
  cov_v <- t(res$y_mat) %*% res$y_mat
  
  max(abs(sum(cov_u) - sum(diag(cov_u))), abs(sum(cov_v) - sum(diag(cov_v)))) < 1e-6
 })
 
 expect_true(all(bool_vec))
})
