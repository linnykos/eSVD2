context("Test optimization")

test_that("opt_esvd works for gaussian", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = 1,
                       library_size_vec = 1)

  ## Determine initialization
  init <- initialize_esvd(dat, k = k, family = "gaussian", nuisance_param_vec = NA,
                          library_size_vec = 1, config = initialization_options())

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = .gaussian,
                  nuisance_param_vec = 1, library_size_vec = NA,
                  max_iter = 15, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "y", "loss"))))
})

test_that("opt_esvd works for poisson", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- log(tcrossprod(x_mat, y_mat))

  ## Simulate data
  dat <- generate_data(nat_mat, family = "poisson", nuisance_param_vec = NA,
                       library_size_vec = 1)

  ## Determine initialization
  init <- initialize_esvd(dat, k = k, family = "poisson", nuisance_param_vec = NA,
                          library_size_vec = 1, config = initialization_options())

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = .poisson,
                  nuisance_param_vec = 1, library_size_vec = NA,
                  max_iter = 15, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "y", "loss"))))
})

test_that("opt_esvd works for negative binomial", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "neg_binom", nuisance_param_vec = 10,
                       library_size_vec = 1)

  ## Determine initialization
  init <- initialize_esvd(dat, k = k, family = "neg_binom", nuisance_param_vec = 10,
                          library_size_vec = 1, config = initialization_options())

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = .neg_binom,
                  nuisance_param_vec = 10, library_size_vec = 1,
                  max_iter = 15, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "y", "loss"))))
  expect_true(all(sign(init$x_mat %*% t(init$y_mat)) == sign(res$x %*% t(res$y))))
})

test_that("opt_esvd works for exponential", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- -tcrossprod(x_mat, y_mat) / 10

  ## Simulate data
  dat <- generate_data(nat_mat, family = "exponential", nuisance_param_vec = NA,
                              library_size_vec = NA)

  ## Determine initialization
  init <- initialize_esvd(dat, k = k, family = "exponential", nuisance_param_vec = NA,
                          library_size_vec = 1, config = initialization_options())

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = .exponential,
                  nuisance_param_vec = NA, library_size_vec = 1,
                  max_iter = 15, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "y", "loss"))))
  expect_true(all(sign(init$x_mat %*% t(init$y_mat)) == sign(res$x %*% t(res$y))))
})

test_that("opt_esvd works for curved gaussian", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2,
                              library_size_vec = 1, tol = 1e-3)

  ## Determine initialization
  init <- initialize_esvd(dat, k = k, family = "curved_gaussian", nuisance_param_vec = 2,
                          library_size_vec = 1, config = initialization_options())

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = .curved_gaussian,
                  nuisance_param_vec = 2, library_size_vec = 1,
                  max_iter = 10, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "y", "loss"))))
  expect_true(all(sign(init$x_mat %*% t(init$y_mat)) == sign(res$x %*% t(res$y))))
})

test_that("opt_esvd works for bernoulli", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "bernoulli")

  ## Determine initialization
  init <- initialize_esvd(dat, k = k, family = "bernoulli",
                          library_size_vec = 1, config = initialization_options())

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = .bernoulli,
                  library_size_vec = 1, max_iter = 10, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x", "y", "loss"))))
})
