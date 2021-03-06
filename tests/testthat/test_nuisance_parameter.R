context("Test nuisance parameter")

test_that("initialize_nuisance_param works for neg_binom", {
  set.seed(5)
  nat_mat <- matrix(-1/5, 100, 5)
  dat <- generate_data(nat_mat, family = "neg_binom", nuisance_param_vec = 10)

  res <- initialize_esvd(dat, k = 1, family = "poisson", library_size_vec = 1)
  res2 <- opt_esvd(res$x_mat, res$y_mat, dat, family = "poisson",
                   library_size_vec = 1)

  nuisance_vec <- initialize_nuisance_param(dat, res2$x %*% t(res2$y), family = "neg_binom",
                                            library_size_vec = res$library_size_vec)

  expect_true(all(nuisance_vec > 0))
  expect_true(is.numeric(nuisance_vec))
  expect_true(length(nuisance_vec) == ncol(dat))
})

test_that("initialize_nuisance_param works for neg_binom with library and nuisance", {
  set.seed(5)
  n <- 500; p <- 5
  nat_mat <- matrix(-1/5, n, p)
  dat <- generate_data(nat_mat, family = "neg_binom",
                       nuisance_param_vec = c(1:p)*10, library_size_vec = 1:n)

  res <- initialize_esvd(dat, k = 1, family = "poisson", library_size_vec = NA)
  res2 <- opt_esvd(res$x_mat, res$y_mat, dat, family = "poisson",
                   library_size_vec = 1:n)

  nuisance_vec <- initialize_nuisance_param(dat, res2$x %*% t(res2$y), family = "neg_binom",
                                            library_size_vec = res2$library_size_vec)

  expect_true(all(nuisance_vec > 0))
  expect_true(is.numeric(nuisance_vec))
  expect_true(length(nuisance_vec) == ncol(dat))
})

test_that("initialize_nuisance_param works for curved_gaussian", {
  set.seed(5)
  nat_mat <- matrix(1/5, 100, 5)
  dat <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2)

  res <- initialize_esvd(dat, k = 1, family = "exponential", library_size_vec = 1)
  res2 <- opt_esvd(res$x_mat, res$y_mat, dat, family = "exponential",
                   library_size_vec = 1)

  nuisance_vec <- initialize_nuisance_param(dat, res2$x %*% t(res2$y), family = "curved_gaussian",
                                            library_size_vec = res$library_size_vec)

  expect_true(all(nuisance_vec > 0))
  expect_true(is.numeric(nuisance_vec))
  expect_true(length(nuisance_vec) == ncol(dat))
})

test_that("initialize_nuisance_param works for curved_gaussian with library and nuisance", {
  set.seed(5)
  n <- 500; p <- 5
  nat_mat <- matrix(1/5, n, p)
  dat <- generate_data(nat_mat, family = "curved_gaussian",
                       nuisance_param_vec = c(1:p)/5, library_size_vec = 1:n)

  res <- initialize_esvd(dat, k = 1, family = "exponential", library_size_vec = NA)
  res2 <- opt_esvd(res$x_mat, res$y_mat, dat, family = "exponential",
                   library_size_vec = res$library_size_vec)

  nuisance_vec <- initialize_nuisance_param(dat, res2$x %*% t(res2$y), family = "curved_gaussian",
                                            library_size_vec = res$library_size_vec)

  expect_true(all(nuisance_vec > 0))
  expect_true(is.numeric(nuisance_vec))
  expect_true(length(nuisance_vec) == ncol(dat))
})

