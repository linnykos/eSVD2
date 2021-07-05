context("Test initialization")

## .matrix_completion is correct

test_that(".matrix_completion works", {
 set.seed(10)
 dat <- abs(matrix(rnorm(40), nrow = 10, ncol = 4))
 dat[sample(1:prod(dim(dat)), 10)] <- NA
 res <- .matrix_completion(dat, k = 2)

 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(dat)))
})

########################

## .determine_initial_matrix is correct

test_that(".determine_initial_matrix works", {
 set.seed(10)
 dat <- matrix(1:40, nrow = 10, ncol = 4)
 family <- "poisson"
 family <- .string_to_distr_funcs(family)

 res <- .determine_initial_matrix(dat, k = 2, family, nuisance_param_vec = NA)

 expect_true(all(sort(names(res)) == sort(c("nat_mat", "domain"))))
 expect_true(length(res$domain) == 2)
 expect_true(res$domain[1] <= res$domain[2])
 expect_true(all(dim(res$nat_mat) == dim(dat)))
})


test_that(".determine_initial_matrix respects max_val", {
 set.seed(10)
 dat <- matrix(1:40, nrow = 10, ncol = 4)
 family <- "poisson"
 family <- .string_to_distr_funcs(family)

 max_val <- 2
 res <- .determine_initial_matrix(dat, k = 2, family, nuisance_param_vec = NA,
                                  max_val = max_val)

 expect_true(all(res$nat_mat >= 0))
 expect_true(all(res$nat_mat <= max_val))
})

################################

## .fix_rank_defficiency is correct

test_that(".fix_rank_defficiency works", {
 set.seed(10)
 x_mat <- abs(matrix(stats::rnorm(12), nrow = 4, ncol = 3))
 y_mat <- abs(matrix(stats::rnorm(15), nrow = 5, ncol = 3))

 res <- .fix_rank_defficiency(x_mat, y_mat, domain = c(0.1, 100))

 expect_true(length(res) == 2)
 expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat"))))

 nat_mat <- res$x_mat %*% t(res$y_mat)
 expect_true(all(nat_mat > 0))
 expect_true(Matrix::rankMatrix(nat_mat) == 3)
})

test_that(".fix_rank_defficiency can actually fix the rank", {
 set.seed(10)
 x_mat <- abs(matrix(stats::rnorm(12), nrow = 4, ncol = 3))
 y_mat <- abs(matrix(stats::rnorm(15), nrow = 5, ncol = 3))
 y_mat[,3] <- 0

 nat_mat <- x_mat %*% t(y_mat)

 res <- .fix_rank_defficiency(x_mat, y_mat, domain = c(min(nat_mat)/2, 100))

 nat_mat <- res$x_mat %*% t(res$y_mat)
 expect_true(all(nat_mat > 0))
 expect_true(Matrix::rankMatrix(nat_mat) == 3)
})

#############################

## initialize_esvd is correct

test_that("initialize_esvd works", {
 set.seed(10)
 dat <- matrix(1:40, nrow = 10, ncol = 4)

 res <- initialize_esvd(dat, k = 2, family = "poisson", library_size_vec = 1)

 expect_true(is.list(res))
 expect_true(class(res) == "eSVD")
 expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat", "domain",
                                            "library_size_vec",
                                            "nuisance_param_vec", "b_mat"))))
 expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
 expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
})

# test_that("initialize_esvd works for gaussian with library size", {
#   set.seed(5)
#   canon_vec <- seq(10, 50, length.out = 5)
#   nat_mat <- sapply(canon_vec, function(x){rep(x, 1000)})
#   dat <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = 1)
#
#   res <- initialize_esvd(dat, k = 2, family = "gaussian", library_size_vec = NA)
#
#   expect_true(class(res) == "eSVD")
#   expect_true(length(res$library_size_vec) == nrow(dat))
#   expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
#   expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
# })

test_that("initialize_esvd works for poisson with library size", {
  set.seed(5)
  canon_vec <- seq(10, 50, length.out = 5)
  nat_mat <- sapply(log(canon_vec), function(x){rep(x, 5000)})
  dat <- generate_data(nat_mat, family = "poisson")

  res <- initialize_esvd(dat, k = 2, family = "poisson", library_size_vec = NA)

  expect_true(class(res) == "eSVD")
  expect_true(length(res$library_size_vec) == nrow(dat))
  expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
  expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
})

# test_that("initialize_esvd works for curved_gaussian with library size", {
#   set.seed(5)
#   canon_vec <- seq(10, 50, length.out = 5)
#   nat_mat <- sapply(1/canon_vec, function(x){rep(x, 1000)})
#   dat <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2)
#
#   res <- initialize_esvd(dat, k = 2, family = "curved_gaussian", library_size_vec = NA,
#                          nuisance_param_vec = 2)
#
#   expect_true(class(res) == "eSVD")
#   expect_true(length(res$library_size_vec) == nrow(dat))
#   expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
#   expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
# })

# test_that("initialize_esvd works for exponential with library size", {
#   set.seed(5)
#   canon_vec <- seq(10, 50, length.out = 5)
#   nat_mat <- sapply(-1/canon_vec, function(x){rep(x, 2000)})
#   dat <- generate_data(nat_mat, family = "exponential")
#
#   res <- initialize_esvd(dat, k = 2, family = "exponential", library_size_vec = NA,
#                          nuisance_param_vec = NA)
#
#   expect_true(class(res) == "eSVD")
#   expect_true(length(res$library_size_vec) == nrow(dat))
#   expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
#   expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
# })

test_that("initialize_esvd works for negative binomial with library size", {
  set.seed(5)
  nat_mat <- matrix(-1/c(1:30), 5, 6)
  dat <- generate_data(nat_mat, family = "neg_binom", nuisance_param_vec = 10)

  res <- initialize_esvd(dat, k = 2, family = "neg_binom", library_size_vec = NA,
                         nuisance_param_vec = 10)

  expect_true(class(res) == "eSVD")
  expect_true(length(res$library_size_vec) == nrow(dat))
  expect_true(all(dim(res$x_mat) == c(nrow(dat), 2)))
  expect_true(all(dim(res$y_mat) == c(ncol(dat), 2)))
})

test_that("initialize_esvd does not suffer from a strange numeric issue with domain", {
 set.seed(10)

 n <- 100; p <- 150; k <- 5
 x_mat <- matrix(abs(stats::rnorm(n*k)), nrow = n, ncol = k)
 y_mat <- matrix(abs(stats::rnorm(p*k)), nrow = p, ncol = k)
 nat_mat <- (x_mat %*% t(y_mat))/10

 dat <- generate_data(nat_mat, family = "poisson")

 init_res <- initialize_esvd(dat, k = k, family = "poisson", nuisance_param_vec = NA, library_size_vec = 1,
                                    config = eSVD2::initialization_options())

 nat_mat <- init_res$x_mat %*% t(init_res$y_mat)
 expect_true(all(nat_mat >= init_res$domain[1]))
 expect_true(all(nat_mat <= init_res$domain[2]))

 ## this is what happens when you don't use .check_domain
 # expect_true(min(nat_mat) >= init_res$domain[1])
})

# test_that("initialize_esvd domain is not check for gaussian", {
#   set.seed(10)
#   dat <- matrix(stats::rnorm(1:40/10), nrow = 10, ncol = 4)
#   expect_true(min(dat) < 0)
#
#   res <- initialize_esvd(dat, k = 2, family = "gaussian", library_size_vec = 1)
#   expect_true(class(res) == "eSVD")
#
#   expect_error(initialize_esvd(dat, k = 2, family = "curved_gaussian"))
# })

test_that("initialize_svd works for missing values", {
  set.seed(123)
  n <- 10; p <- 15; k <- 2
  nuisance_param_vec <- runif(p, 0, 5)
  u_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  v_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(u_mat, v_mat)
  dat <- eSVD2::generate_data(
    nat_mat, family = "poisson", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = 1
  )
  dat[sample(1:prod(dim(dat)), 5)] <- NA

  res <- initialize_esvd(dat, k = k, family = "poisson",
                         nuisance_param_vec = nuisance_param_vec,
                         library_size_vec = 1)

  expect_true(is.list(res))
})
