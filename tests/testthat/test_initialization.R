context("Test initialization")

## initialize_esvd is correct

test_that("initialize_esvd works", {
 set.seed(10)
 dat <- matrix(1:40, nrow = 10, ncol = 4)

 res <- initialize_esvd(dat, k = 2, family = "poisson")

 expect_true(is.list(res))
 expect_true(class(res) == "eSVD")
 expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat", "covariates",
                                            "offset_vec",
                                            "nuisance_param_vec", "b_mat"))))
 expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
 expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
})

test_that("initialize_esvd works for poisson with covariates", {
  set.seed(5)
  canon_vec <- 1:50
  nat_mat <- sapply(log(canon_vec)+1, function(x){rep(x, 100)})
  dat <- generate_data(nat_mat, family = "poisson")

  res <- initialize_esvd(dat, k = 2, family = "poisson",
                         covariates = matrix(1, nrow = nrow(dat), ncol = 1))

  expect_true(class(res) == "eSVD")
  expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
  expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
})

test_that("initialize_esvd works for negative binomial", {
  set.seed(5)
  nat_mat <- matrix(1/c(1:30), 5, 6)
  dat <- generate_data(nat_mat, family = "neg_binom2", nuisance_param_vec = 10)

  res <- initialize_esvd(dat, k = 2, family = "neg_binom2")

  expect_true(class(res) == "eSVD")
  expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
  expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
})

test_that("initialize_esvd works for missing values", {
  set.seed(123)
  n <- 10; p <- 15; k <- 2
  nuisance_param_vec <- runif(p, 0, 5)
  u_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  v_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(u_mat, v_mat)
  dat <- eSVD2::generate_data(
    nat_mat, family = "poisson"
  )
  dat[sample(1:prod(dim(dat)), 5)] <- NA

  res <- initialize_esvd(dat, k = k, family = "poisson")

  expect_true(is.list(res))
})
