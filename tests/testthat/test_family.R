context("Test distribution family functions")

run_test <- function(
  dat, x_mat, y_mat, family, offset = rnorm(nrow(x_mat)),
  nuisance_param_vec = rep(1, nrow(y_mat)), library_size_vec = rep(1, nrow(x_mat)),
  l2pen = list(x = 0.001, y = 0.002, z = 0.003), tol = 1e-6
)
{
  n <- nrow(dat)
  p <- ncol(dat)
  k <- ncol(x_mat)

  c_mat <- matrix(offset, ncol = 1)
  z_mat <- matrix(1, nrow = p, ncol = 1)
  xc_mat <- cbind(x_mat, c_mat)
  yz_mat <- cbind(y_mat, z_mat)
  loader <- data_loader(dat)
  family <- esvd_family(family)
  library_size_vec <- .parse_library_size(dat, library_size_vec)
  if(length(nuisance_param_vec) == 1)
    nuisance_param_vec <- rep(nuisance_param_vec, p)
  yzind <- (1:k) - 1

  ############# Test objective functions #############
  loss1 <- objfn_all_r(
    XC = xc_mat, YZ = yz_mat, k = k, loader = loader, family = family,
    s = library_size_vec, gamma = nuisance_param_vec,
    l2penx = l2pen$x, l2peny = l2pen$y, l2penz = l2pen$z
  )

  loss2 <- numeric(n)
  for(i in 1:n)
  {
    loss2[i] <- objfn_Xi_r(
      XCi = xc_mat[i, ], YZ = yz_mat, k = k, loader = loader, row_ind = i - 1,
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec,
      l2penx = l2pen$x
    )
    loss2[i] <- loss2[i] * sum(!is.na(dat[i, ]))
  }

  loss3 <- numeric(p)
  for(j in 1:p)
  {
    loss3[j] <- objfn_YZj_r(
      XC = xc_mat, YZj = yz_mat[j, ], k = k, YZind = yzind,
      loader = loader, col_ind = j - 1, family = family,
      s = library_size_vec, gammaj = nuisance_param_vec[j],
      l2peny = l2pen$y, l2penz = l2pen$z
    )
    loss3[j] <- loss3[j] * sum(!is.na(dat[, j]))
  }

  # The three loss values should be equal
  loss2 <- sum(loss2) + l2pen$y * sum(y_mat^2) + l2pen$z * sum(z_mat^2)
  loss2 <- loss2 / sum(!is.na(dat))
  loss3 <- sum(loss3) + l2pen$x * sum(x_mat^2)
  loss3 <- loss3 / sum(!is.na(dat))

  expect_lt(abs(loss1 - loss2), tol)
  expect_lt(abs(loss1 - loss3), tol)

  ############# Test gradients #############
  for(i in 1:n)
  {
    grad1 <- grad_Xi_r(
      XCi = xc_mat[i, ], YZ = yz_mat, k = k, loader = loader, row_ind = i - 1,
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec,
      l2penx = l2pen$x
    )
    grad2 <- numDeriv::grad(
      function(x, ...) objfn_Xi_r(XCi = c(x, c_mat[i, ]), ...),
      x_mat[i, ], side = NULL, method.args = list(r = 6),
      YZ = yz_mat, k = k, loader = loader, row_ind = i - 1,
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec,
      l2penx = l2pen$x
    )
    expect_lt(max(abs(grad1 - grad2)), tol)
  }

  for(j in 1:p)
  {
    grad1 <- grad_YZj_r(
      XC = xc_mat, YZj = yz_mat[j, ], k = k, YZind = yzind,
      loader = loader, col_ind = j - 1, family = family,
      s = library_size_vec, gammaj = nuisance_param_vec[j],
      l2peny = l2pen$y, l2penz = l2pen$z
    )
    grad2 <- numDeriv::grad(
      function(x, ...) { yzj <- yz_mat[j, ]; yzj[yzind + 1] <- x; objfn_YZj_r(YZj = yzj, ...) },
      yz_mat[j, yzind + 1], side = NULL, method.args = list(r = 6),
      XC = xc_mat, , k = k, YZind = yzind,
      loader = loader, col_ind = j - 1, family = family,
      s = library_size_vec, gammaj = nuisance_param_vec[j],
      l2peny = l2pen$y, l2penz = l2pen$z
    )
    expect_lt(max(abs(grad1 - grad2)), tol)
  }

  ############# Test Hessians #############
  for(i in 1:n)
  {
    hess1 <- hessian_Xi_r(
      XCi = xc_mat[i, ], YZ = yz_mat, k = k, loader = loader, row_ind = i - 1,
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec,
      l2penx = l2pen$x
    )
    hess2 <- numDeriv::hessian(
      function(x, ...) objfn_Xi_r(XCi = c(x, c_mat[i, ]), ...),
      x_mat[i, ],
      YZ = yz_mat, k = k, loader = loader, row_ind = i - 1,
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec,
      l2penx = l2pen$x
    )
    expect_lt(max(abs(hess1 - hess2)), tol)
  }

  for(j in 1:p)
  {
    hess1 <- hessian_YZj_r(
      XC = xc_mat, YZj = yz_mat[j, ], k = k, YZind = yzind,
      loader = loader, col_ind = j - 1, family = family,
      s = library_size_vec, gammaj = nuisance_param_vec[j],
      l2peny = l2pen$y, l2penz = l2pen$z
    )
    hess2 <- numDeriv::hessian(
      function(x, ...) { yzj <- yz_mat[j, ]; yzj[yzind + 1] <- x; objfn_YZj_r(YZj = yzj, ...) },
      yz_mat[j, yzind + 1],
      XC = xc_mat, , k = k, YZind = yzind,
      loader = loader, col_ind = j - 1, family = family,
      s = library_size_vec, gammaj = nuisance_param_vec[j],
      l2peny = l2pen$y, l2penz = l2pen$z
    )
    expect_lt(max(abs(hess1 - hess2)), tol)
  }
}


# ######################## Gaussian ########################

test_that("Functions for Gaussian distribution", {
  # Simulate natural parameter matrix
  set.seed(123)
  n <- 10
  p <- 15
  k <- 2
  nuisance_param_vec <- runif(p, 0, 5)
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Simulate data with default library size (all one)
  dat <- eSVD2::generate_data(
    nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = 1
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1)

  # Simulate data with a library size vector
  library_size_vec <- sample(10:20, n, replace = TRUE)
  dat <- eSVD2::generate_data(
    nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = library_size_vec
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec)
})

# ######################## Curved Gaussian ########################

test_that("Functions for curved-Gaussian distribution", {
  # Simulate natural parameter matrix
  set.seed(123)
  n <- 10
  p <- 15
  k <- 2
  nuisance_param_vec <- runif(p, 0, 5)
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Simulate data with default library size (all one)
  dat <- eSVD2::generate_data(
    nat_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = 1, tol = 1e-3
  )

  # Force a positive offset
  offset = runif(n)

  # Test
  run_test(dat, x_mat, y_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1, offset = offset, tol = 1e-5)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1, offset = offset, tol = 1e-5)

  # Simulate data with a library size vector
  library_size_vec <- sample(1:10, n, replace = TRUE)
  dat <- eSVD2::generate_data(
    nat_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = library_size_vec, tol = 1e-3
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec, offset = offset, tol = 1e-4)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec, offset = offset, tol = 1e-4)
})

# ######################## Exponential ########################

test_that("Functions for exponential distribution", {
  # Simulate natural parameter matrix
  set.seed(123)
  n <- 10
  p <- 15
  k <- 2
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Simulate data with default library size (all one)
  dat <- eSVD2::generate_data(
    nat_mat, family = "exponential", nuisance_param_vec = NA,
    library_size_vec = 1
  )

  # Force a negative offset
  offset = -runif(n)

  # Test
  run_test(dat, x_mat, y_mat, family = "exponential", nuisance_param_vec = NA,
           library_size_vec = 1, offset = offset)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "exponential", nuisance_param_vec = NA,
           library_size_vec = 1, offset = offset)

  # Simulate data with a library size vector
  library_size_vec <- sample(10:20, n, replace = TRUE)
  dat <- eSVD2::generate_data(
    nat_mat, family = "exponential", nuisance_param_vec = NA,
    library_size_vec = library_size_vec
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "exponential", nuisance_param_vec = NA,
           library_size_vec = library_size_vec, offset = offset)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "exponential", nuisance_param_vec = NA,
           library_size_vec = library_size_vec, offset = offset)
})

######################## Poisson ########################

test_that("Functions for Poisson distribution", {
  # Simulate natural parameter matrix
  set.seed(123)
  n <- 10
  p <- 15
  k <- 2
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  # Simulate data with default library size (all one)
  dat <- eSVD2::generate_data(
    nat_mat, family = "poisson", nuisance_param_vec = NA,
    library_size_vec = 1
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "poisson", nuisance_param_vec = NA,
           library_size_vec = 1)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "poisson", nuisance_param_vec = NA,
           library_size_vec = 1)

  # Simulate data with a library size vector
  library_size_vec <- sample(10:20, n, replace = TRUE)
  dat <- eSVD2::generate_data(
    nat_mat, family = "poisson", nuisance_param_vec = NA,
    library_size_vec = library_size_vec
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "poisson", nuisance_param_vec = NA,
           library_size_vec = library_size_vec)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "poisson", nuisance_param_vec = NA,
           library_size_vec = library_size_vec)
})

######################## Negative binomial ########################

test_that("Functions for negative binomial distribution", {
  # Simulate data
  set.seed(123)
  n <- 10
  p <- 15
  k <- 2
  nuisance_param_vec <- runif(p, 0, 10)
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  dat <- eSVD2::generate_data(
    nat_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = 1
  )

  # Force a negative offset
  offset = -runif(n)

  # Test
  run_test(dat, x_mat, y_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1, offset = offset)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1, offset = offset)

  # Simulate data with a library size vector
  library_size_vec <- sample(10:20, n, replace = TRUE)
  dat <- eSVD2::generate_data(
    nat_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = library_size_vec
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec, offset = offset, tol = 1e-5)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec, offset = offset, tol = 1e-5)
})

######################## Negative binomial 2 ########################

test_that("Functions for negative binomial 2 distribution", {
  # Simulate data
  set.seed(123)
  n <- 10
  p <- 15
  k <- 2
  nuisance_param_vec <- runif(p, 0, 10)
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  dat <- eSVD2::generate_data(
    nat_mat, family = "neg_binom2", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = 1
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "neg_binom2", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "neg_binom2", nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1)
})

# ######################## Bernoulli ########################

test_that("Functions for Bernoulli distribution", {
  # Simulate data
  set.seed(123)
  n <- 10
  p <- 15
  k <- 2
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  dat <- eSVD2::generate_data(
    nat_mat, family = "bernoulli", nuisance_param_vec = NA,
    library_size_vec = 1
  )

  # Test
  run_test(dat, x_mat, y_mat, family = "bernoulli", nuisance_param_vec = NA,
           library_size_vec = 1)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, family = "bernoulli", nuisance_param_vec = NA,
           library_size_vec = 1)
})
