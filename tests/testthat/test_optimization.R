context("Test optimization")

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
  init <- initialize_esvd(dat, k = k, family = "poisson")

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = "poisson",
                  nuisance_param_vec = 1, library_size_vec = 1,
                  max_iter = 15, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat", "b_mat", "loss",
                                             "nuisance_param_vec", "library_size_vec",
                                             "offset_vec", "covariates"))))
})

test_that("opt_esvd works for negative binomial2", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  ## Simulate data
  dat <- generate_data(nat_mat, family = "neg_binom2", nuisance_param_vec = 10,
                       library_size_vec = 1)

  ## Determine initialization
  init <- initialize_esvd(dat, k = k, family = "neg_binom2")

  ## Optimization
  res <- opt_esvd(init$x_mat, init$y_mat, dat, family = "neg_binom2",
                  nuisance_param_vec = 10, library_size_vec = 1,
                  max_iter = 15, verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat", "b_mat", "loss",
                                             "nuisance_param_vec", "library_size_vec",
                                             "offset_vec", "covariates"))))
})
