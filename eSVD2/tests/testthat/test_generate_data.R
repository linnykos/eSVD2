context("Testing generate data")

## generate_data is correct

test_that("generate_data works", {
 set.seed(10)
 nat_mat <- matrix(log(c(11:40)), 5, 6)
 res <- generate_data(nat_mat, family = "poisson")
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data works for poisson", {
 set.seed(5)
 nat_mat <- matrix(log(c(11:40)), 5, 6)
 res <- generate_data(nat_mat, family = "poisson")
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data works for curved gaussian", {
 set.seed(5)
 nat_mat <- matrix(1/c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data works for exponential", {
 set.seed(5)
 nat_mat <- matrix(-c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "exponential")
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data works for neg_binom", {
 set.seed(5)
 nat_mat <- matrix(-1/c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "neg_binom", nuisance_param_vec = 10)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})