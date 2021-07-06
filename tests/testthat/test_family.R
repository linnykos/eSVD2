context("Test distribution family functions")

run_test <- function(dat, x_mat, y_mat, family, nuisance_param_vec = rep(1, nrow(y_mat)), library_size_vec = rep(1, nrow(x_mat)))
{
  n <- nrow(dat)
  p <- ncol(dat)
  library_size_vec <- .parse_library_size(dat, library_size_vec)

  # Test objective functions
  loss1 <- objfn_all(
    X = x_mat, Y = y_mat, B = NULL, Z = NULL, A = dat,
    family = family, s = library_size_vec, gamma = nuisance_param_vec
  )
  loss2 <- numeric(n)
  for(i in 1:n)
  {
    loss2[i] <- objfn_Xi(
      Xi = x_mat[i, ], Y = y_mat, B = NULL, Zi = NULL, Ai = dat[i, ],
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec
    )
    loss2[i] <- loss2[i] * sum(!is.na(dat[i, ]))
  }
  loss3 <- numeric(p)
  for(j in 1:p)
  {
    loss3[j] <- objfn_Yj(
      Yj = y_mat[j, ], X = x_mat, Bj = NULL, Z = NULL, Aj = dat[, j],
      family = family, s = library_size_vec, gammaj = nuisance_param_vec[j]
    )
    loss3[j] <- loss3[j] * sum(!is.na(dat[, j]))
  }
  # The three loss values should be equal
  expect_lt(loss1 - sum(loss2) / sum(!is.na(dat)), 1e-8)
  expect_lt(loss1 - sum(loss3) / sum(!is.na(dat)), 1e-8)

  # Test gradients
  for(i in 1:n)
  {
    grad1 <- grad_Xi(
      Xi = x_mat[i, ], Y = y_mat, B = NULL, Zi = NULL, Ai = dat[i, ],
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec
    )
    grad2 <- numDeriv::grad(
      objfn_Xi, x_mat[i, ], side = NULL, Y = y_mat, B = NULL, Zi = NULL, Ai = dat[i, ],
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec
    )
    expect_lt(max(abs(grad1 - grad2)), 1e-6)
  }
  for(j in 1:p)
  {
    grad1 <- grad_Yj(
      Yj = y_mat[j, ], X = x_mat, Bj = NULL, Z = NULL, Aj = dat[, j],
      family = family, s = library_size_vec, gammaj = nuisance_param_vec[j]
    )
    grad2 <- numDeriv::grad(
      objfn_Yj, y_mat[j, ], side = NULL, X = x_mat, Bj = NULL, Z = NULL, Aj = dat[, j],
      family = family, s = library_size_vec, gammaj = nuisance_param_vec[j]
    )
    expect_lt(max(abs(grad1 - grad2)), 1e-6)
  }

  # Test Hessians
  for(i in 1:n)
  {
    hess1 <- hessian_Xi(
      Xi = x_mat[i, ], Y = y_mat, B = NULL, Zi = NULL, Ai = dat[i, ],
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec
    )
    hess2 <- numDeriv::hessian(
      objfn_Xi, x_mat[i, ], Y = y_mat, B = NULL, Zi = NULL, Ai = dat[i, ],
      family = family, si = library_size_vec[i], gamma = nuisance_param_vec
    )
    expect_lt(max(abs(hess1 - hess2)), 1e-5)
  }
  for(j in 1:p)
  {
    hess1 <- hessian_Yj(
      Yj = y_mat[j, ], X = x_mat, Bj = NULL, Z = NULL, Aj = dat[, j],
      family = family, s = library_size_vec, gammaj = nuisance_param_vec[j]
    )
    hess2 <- numDeriv::hessian(
      objfn_Yj, y_mat[j, ], X = x_mat, Bj = NULL, Z = NULL, Aj = dat[, j],
      family = family, s = library_size_vec, gammaj = nuisance_param_vec[j]
    )
    expect_lt(max(abs(hess1 - hess2)), 1e-5)
  }
}


# ######################## Gaussian ########################
#
# test_that("Functions for Gaussian distribution", {
#   # Simulate natural parameter matrix
#   set.seed(123)
#   n <- 10
#   p <- 15
#   k <- 2
#   nuisance_param_vec <- runif(p, 0, 5)
#   x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
#   y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
#   nat_mat <- tcrossprod(x_mat, y_mat)
#
#   # Simulate data with default library size (all one)
#   dat <- eSVD2::generate_data(
#     nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
#     library_size_vec = 1
#   )
#
#   # Test
#   run_test(dat, x_mat, y_mat, .esvd.gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = 1)
#
#   # Test missing values
#   dat[sample(length(dat), n * p * 0.1)] <- NA
#   run_test(dat, x_mat, y_mat, .esvd.gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = 1)
#
#   # Simulate data with a library size vector
#   library_size_vec <- sample(10:20, n, replace = TRUE)
#   dat <- eSVD2::generate_data(
#     nat_mat, family = "gaussian", nuisance_param_vec = nuisance_param_vec,
#     library_size_vec = library_size_vec
#   )
#
#   # Test
#   run_test(dat, x_mat, y_mat, .esvd.gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = library_size_vec)
#
#   # Test missing values
#   dat[sample(length(dat), n * p * 0.1)] <- NA
#   run_test(dat, x_mat, y_mat, .esvd.gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = library_size_vec)
# })
#
# ######################## Curved Gaussian ########################
#
# test_that("Functions for curved-Gaussian distribution", {
#   # Simulate natural parameter matrix
#   set.seed(123)
#   n <- 10
#   p <- 15
#   k <- 2
#   nuisance_param_vec <- runif(p, 0, 5)
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- tcrossprod(x_mat, y_mat)
#
#   # Simulate data with default library size (all one)
#   dat <- eSVD2::generate_data(
#     nat_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
#     library_size_vec = 1, tol = 1e-3
#   )
#
#   # Test
#   run_test(dat, x_mat, y_mat, .esvd.curved_gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = 1)
#
#   # Test missing values
#   dat[sample(length(dat), n * p * 0.1)] <- NA
#   run_test(dat, x_mat, y_mat, .esvd.curved_gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = 1)
#
#   # Simulate data with a library size vector
#   library_size_vec <- sample(10:20, n, replace = TRUE)
#   dat <- eSVD2::generate_data(
#     nat_mat, family = "curved_gaussian", nuisance_param_vec = nuisance_param_vec,
#     library_size_vec = library_size_vec, tol = 1e-3
#   )
#
#   # Test
#   run_test(dat, x_mat, y_mat, .esvd.curved_gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = library_size_vec)
#
#   # Test missing values
#   dat[sample(length(dat), n * p * 0.1)] <- NA
#   run_test(dat, x_mat, y_mat, .esvd.curved_gaussian, nuisance_param_vec = nuisance_param_vec,
#            library_size_vec = library_size_vec)
# })
#
# ######################## Exponential ########################
#
# test_that("Functions for exponential distribution", {
#   # Simulate natural parameter matrix
#   set.seed(123)
#   n <- 10
#   p <- 15
#   k <- 2
#   x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
#   y_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
#   nat_mat <- tcrossprod(x_mat, y_mat)
#
#   # Simulate data with default library size (all one)
#   dat <- eSVD2::generate_data(
#     nat_mat, family = "exponential", nuisance_param_vec = NA,
#     library_size_vec = 1
#   )
#
#   # Test
#   run_test(dat, x_mat, y_mat, .esvd.exponential, nuisance_param_vec = NA,
#            library_size_vec = 1)
#
#   # Test missing values
#   dat[sample(length(dat), n * p * 0.1)] <- NA
#   run_test(dat, x_mat, y_mat, .esvd.exponential, nuisance_param_vec = NA,
#            library_size_vec = 1)
#
#   # Simulate data with a library size vector
#   library_size_vec <- sample(10:20, n, replace = TRUE)
#   dat <- eSVD2::generate_data(
#     nat_mat, family = "exponential", nuisance_param_vec = NA,
#     library_size_vec = library_size_vec
#   )
#
#   # Test
#   run_test(dat, x_mat, y_mat, .esvd.exponential, nuisance_param_vec = NA,
#            library_size_vec = library_size_vec)
#
#   # Test missing values
#   dat[sample(length(dat), n * p * 0.1)] <- NA
#   run_test(dat, x_mat, y_mat, .esvd.exponential, nuisance_param_vec = NA,
#            library_size_vec = library_size_vec)
# })

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
  run_test(dat, x_mat, y_mat, .esvd.poisson, nuisance_param_vec = NA,
           library_size_vec = 1)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, .esvd.poisson, nuisance_param_vec = NA,
           library_size_vec = 1)

  # Simulate data with a library size vector
  library_size_vec <- sample(10:20, n, replace = TRUE)
  dat <- eSVD2::generate_data(
    nat_mat, family = "poisson", nuisance_param_vec = NA,
    library_size_vec = library_size_vec
  )

  # Test
  run_test(dat, x_mat, y_mat, .esvd.poisson, nuisance_param_vec = NA,
           library_size_vec = library_size_vec)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, .esvd.poisson, nuisance_param_vec = NA,
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

  # Test
  run_test(dat, x_mat, y_mat, .esvd.neg_binom, nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, .esvd.neg_binom, nuisance_param_vec = nuisance_param_vec,
           library_size_vec = 1)

  # Simulate data with a library size vector
  library_size_vec <- sample(10:20, n, replace = TRUE)
  dat <- eSVD2::generate_data(
    nat_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = library_size_vec
  )

  # Test
  run_test(dat, x_mat, y_mat, .esvd.neg_binom, nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec)

  # Test missing values
  dat[sample(length(dat), n * p * 0.1)] <- NA
  run_test(dat, x_mat, y_mat, .esvd.neg_binom, nuisance_param_vec = nuisance_param_vec,
           library_size_vec = library_size_vec)
})

# ######################## Bernoulli ########################
#
# test_that("Functions for Bernoulli distribution", {
#   # Simulate data
#   set.seed(123)
#   n <- 10
#   p <- 15
#   k <- 2
#   x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
#   y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
#   nat_mat <- tcrossprod(x_mat, y_mat)
#   dat <- eSVD2::generate_data(
#     nat_mat, family = "bernoulli", nuisance_param_vec = NA,
#     library_size_vec = 1
#   )
#
#   # Test
#   run_test(dat, x_mat, y_mat, .esvd.bernoulli, nuisance_param_vec = NA,
#            library_size_vec = 1)
#
#   # Test missing values
#   dat[sample(length(dat), n * p * 0.1)] <- NA
#   run_test(dat, x_mat, y_mat, .esvd.bernoulli, nuisance_param_vec = NA,
#            library_size_vec = 1)
# })
