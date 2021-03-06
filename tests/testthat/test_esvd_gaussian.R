context("Test Gaussian with fixed variance")

## .evaluate_objective.gaussian is correct

test_that(".evaluate_objective.gaussian works and is correct", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)/(p/2)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
                       library_size_vec = library_size_vec)

  res1 <- .evaluate_objective.gaussian(dat, x_mat, y_mat, nuisance_param_vec = rep(1, p),
                                       library_size_vec = rep(1, n))
  res2 <- .evaluate_objective.gaussian(dat, x_mat, y_mat, nuisance_param_vec = nuisance_param_vec,
                                       library_size_vec = library_size_vec)

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))

  res1b <- 0
  for(i in 1:n) {
    for(j in 1:p) {
      res1b <- res1b + (x_mat[i, ] %*% y_mat[j, ] - dat[i, j])^2 / (n * p)
    }
  }

  res2b <- 0
  for(i in 1:n) {
    for(j in 1:p) {
      res2b <- res2b + (x_mat[i, ] %*% y_mat[j, ] - dat[i, j] / library_size_vec[i])^2 / nuisance_param_vec[j]^2 / (n * p)
    }
  }

  expect_true(abs(res1 - res1b) <= 1e-4)
  expect_true(abs(res2 - res2b) <= 1e-4)
})

#################

## .evaluate_objective_single.gaussian is correct

test_that(".evaluate_objective_single.gaussian works and is correct in the x direction", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)/(p/2)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
                       library_size_vec = library_size_vec)

  i <- 1
  res1 <- .evaluate_objective_single.gaussian(x_mat[i, ], y_mat, dat[i, ], nuisance_param_vec = rep(1, p),
                                              library_size_vec = rep(1, n)[i])
  res2 <- .evaluate_objective_single.gaussian(x_mat[i, ], y_mat, dat[i, ], nuisance_param_vec = nuisance_param_vec,
                                              library_size_vec = library_size_vec[i])

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))

  res1b <- 0
  for(j in 1:p) {
    res1b <- res1b + (x_mat[i, ] %*% y_mat[j, ] - dat[i, j])^2 / p
  }
  res2b <- 0
  for(j in 1:p) {
    res2b <- res2b + (x_mat[i, ] %*% y_mat[j, ] - dat[i, j] / library_size_vec[i])^2 / nuisance_param_vec[j]^2 / p
  }

  expect_true(abs(res1 - res1b) <= 1e-4)
  expect_true(abs(res2 - res2b) <= 1e-4)
})

test_that(".evaluate_objective_single.gaussian works and is correct in the y direction", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)/(p/2)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
                       library_size_vec = library_size_vec)

  j <- 1
  res1 <- .evaluate_objective_single.gaussian(y_mat[j, ], x_mat, dat[, j], nuisance_param_vec = rep(1, p)[j],
                                              library_size_vec = rep(1, n))
  res2 <- .evaluate_objective_single.gaussian(y_mat[j, ], x_mat, dat[, j], nuisance_param_vec = nuisance_param_vec[j],
                                              library_size_vec = library_size_vec)

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))

  res1b <- 0
  for(i in 1:n) {
    res1b <- res1b + (x_mat[i, ] %*% y_mat[j, ] - dat[i, j])^2 / n
  }
  res2b <- 0
  for(i in 1:n) {
    res2b <- res2b + (x_mat[i, ] %*% y_mat[j, ] - dat[i, j] / library_size_vec[i])^2 / nuisance_param_vec[j]^2 / n
  }

  expect_true(abs(res1 - res1b) <= 1e-4)
  expect_true(abs(res2 - res2b) <= 1e-4)
})

###################################

## .gradient_vec.gaussian is correct

test_that(".gradient_vec.gaussian works and is correct in the x direction", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)/(p/2)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
                       library_size_vec = library_size_vec)

  i <- 1
  res1 <- .gradient_vec.gaussian(x_mat[i, ], y_mat, dat[i, ], nuisance_param_vec = rep(1, p),
                                 library_size_vec = rep(1, n)[i])
  res2 <- .gradient_vec.gaussian(x_mat[i, ], y_mat, dat[i, ], nuisance_param_vec = nuisance_param_vec,
                                 library_size_vec = library_size_vec[i])

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))

  res1b <- rep(0, k)
  for(j in 1:p) {
    res1b <- res1b + 2 * c(x_mat[i, ] %*% y_mat[j, ] - dat[i, j]) * y_mat[j, ] / p
  }
  res2b <- rep(0, k)
  for(j in 1:p) {
    res2b <- res2b + 2 * c(x_mat[i, ] %*% y_mat[j, ] - dat[i, j] / library_size_vec[i]) / nuisance_param_vec[j]^2 * y_mat[j,] / p
  }

  expect_true(sum(abs(res1 - res1b)) <= 1e-4)
  expect_true(sum(abs(res2 - res2b)) <= 1e-4)
})

test_that(".gradient_vec.gaussian works and is correct in the y direction", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)/(p/2)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
                       library_size_vec = library_size_vec)

  j <- 1
  res1 <- .gradient_vec.gaussian(y_mat[j, ], x_mat, dat[, j], nuisance_param_vec = rep(1, p)[j],
                                 library_size_vec = rep(1, n))
  res2 <- .gradient_vec.gaussian(y_mat[j, ], x_mat, dat[, j], nuisance_param_vec = nuisance_param_vec[j],
                                 library_size_vec = library_size_vec)

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))

  res1b <- rep(0, k)
  for(i in 1:n) {
    res1b <- res1b + 2 * c(x_mat[i, ] %*% y_mat[j, ] - dat[i, j]) * x_mat[i, ] / n
  }
  res2b <- rep(0, k)
  for(i in 1:n) {
    res2b <- res2b + 2 * c(x_mat[i, ] %*% y_mat[j, ] - dat[i, j] / library_size_vec[i]) / nuisance_param_vec[j]^2 * x_mat[i, ] / n
  }

  expect_true(sum(abs(res1 - res1b)) <= 1e-4)
  expect_true(sum(abs(res2 - res2b)) <= 1e-4)
})

########################################

## .hessian_vec.gaussian is correct

test_that(".hessian_vec.gaussian works and is correct in the x direction", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)/(p/2)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
                       library_size_vec = library_size_vec)

  i <- 1
  res1 <- .hessian_vec.gaussian(x_mat[i, ], y_mat, dat[i, ], nuisance_param_vec = rep(1, p),
                                library_size_vec = rep(1, n)[i])
  res2 <- .hessian_vec.gaussian(x_mat[i, ], y_mat, dat[i, ], nuisance_param_vec = nuisance_param_vec,
                                library_size_vec = library_size_vec[i])

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))

  res1b <- matrix(0, k, k)
  for(j in 1:p) {
    res1b <- res1b + 2 * tcrossprod(y_mat[j, ], y_mat[j, ]) / p
  }
  res2b <- matrix(0, k, k)
  for(j in 1:p) {
    res2b <- res2b + 2 * tcrossprod(y_mat[j, ], y_mat[j, ]) / (p * nuisance_param_vec[j]^2)
  }

  expect_true(sum(abs(res1 - res1b)) <= 1e-4)
  expect_true(sum(abs(res2 - res2b)) <= 1e-4)
})

test_that(".hessian_vec.gaussian works and is correct in the y direction", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)/(p/2)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
                       library_size_vec = library_size_vec)

  j <- 1
  res1 <- .hessian_vec.gaussian(y_mat[j, ], x_mat, dat[, j], nuisance_param_vec = rep(1, p)[j],
                                library_size_vec = rep(1, n))
  res2 <- .hessian_vec.gaussian(y_mat[j, ], x_mat, dat[, j], nuisance_param_vec = nuisance_param_vec[j],
                                library_size_vec = library_size_vec)

  expect_true(is.numeric(res1))
  expect_true(is.numeric(res2))

  res1b <- matrix(0, k, k)
  for(i in 1:n) {
    res1b <- res1b + 2 * tcrossprod(x_mat[i, ], x_mat[i, ]) / n
  }
  res2b <- matrix(0, k, k)
  for(i in 1:n) {
    res2b <- res2b + 2 * tcrossprod(x_mat[i, ], x_mat[i, ]) / (n * nuisance_param_vec[j]^2)
  }

  expect_true(sum(abs(res1 - res1b)) <= 1e-4)
  expect_true(sum(abs(res2 - res2b)) <= 1e-4)
})
