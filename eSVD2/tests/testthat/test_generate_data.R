context("Testing generate data")

## generate_data is correct

test_that("generate_data works", {
 set.seed(10)
 nat_mat <- matrix(log(c(11:40)), 5, 6)
 res <- generate_data(nat_mat, family = "poisson")
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

############

test_that("generate_data works for gaussian", {
 set.seed(5)
 nat_mat <- matrix(log(c(11:40)), 5, 6)
 res <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = 1)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data has the correct mean and variance for gaussian", {
 set.seed(5)
 canon_vec <- seq(10, 50, length.out = 5)
 nat_mat <- sapply(canon_vec, function(x){rep(x, 1000)})
 res <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = 1)
 
 mean_vec <- colMeans(res)
 expect_true(all(abs(mean_vec - canon_vec) <= 2))
 
 var_vec <- apply(res, 2, stats::var)
 expect_true(all(abs(var_vec - 1) <= 2))
})


############

test_that("generate_data works for poisson", {
 set.seed(5)
 nat_mat <- matrix(log(c(11:40)), 5, 6)
 res <- generate_data(nat_mat, family = "poisson")
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data has the correct mean and variance for poisson", {
 set.seed(5)
 canon_vec <- seq(10, 50, length.out = 5)
 nat_mat <- sapply(log(canon_vec), function(x){rep(x, 5000)})
 res <- generate_data(nat_mat, family = "poisson")
 
 mean_vec <- colMeans(res)
 expect_true(all(abs(mean_vec - canon_vec) <= 2))
 
 var_vec <- apply(res, 2, stats::var)
 expect_true(all(abs(var_vec - canon_vec) <= 5))
})

############

test_that("generate_data works for curved gaussian", {
 set.seed(5)
 nat_mat <- matrix(1/c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data has the correct mean and variance for curved guassian", {
 set.seed(5)
 canon_vec <- seq(10, 50, length.out = 5)
 nat_mat <- sapply(1/canon_vec, function(x){rep(x, 1000)})
 res <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2)
 
 mean_vec <- colMeans(res)
 expect_true(all(abs(mean_vec - canon_vec) <= 2))
 
 sd_vec <- apply(res, 2, stats::sd)
 expect_true(all(abs(sd_vec - canon_vec/2) <= 2))
})

############

test_that("generate_data works for exponential", {
 set.seed(5)
 nat_mat <- matrix(-c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "exponential")
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data has the correct mean and variance for curved guassian", {
 set.seed(5)
 canon_vec <- seq(10, 50, length.out = 5)
 nat_mat <- sapply(-1/canon_vec, function(x){rep(x, 2000)})
 res <- generate_data(nat_mat, family = "exponential")
 
 mean_vec <- colMeans(res)
 expect_true(all(abs(mean_vec - canon_vec) <= 2))
 
 sd_vec <- apply(res, 2, stats::sd)
 expect_true(all(abs(sd_vec - canon_vec) <= 3))
})

############

test_that("generate_data works for neg_binom", {
 set.seed(5)
 nat_mat <- matrix(-1/c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "neg_binom", nuisance_param_vec = 10)
 
 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data has the correct mean and variance for neg_binom", {
 set.seed(5)
 canon_vec <- -1/seq(10, 50, length.out = 5)
 nat_mat <- sapply(canon_vec, function(x){rep(x, 2000)})
 nuisance_param_val <- 10
 res <- generate_data(nat_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_val)
 
 mean_vec <- colMeans(res)
 expect_true(all(abs(mean_vec - nuisance_param_val*exp(canon_vec)/(1-exp(canon_vec))) <= 5))
 
 sd_vec <- apply(res, 2, stats::sd)
 expect_true(all(abs(sd_vec - sqrt(nuisance_param_val*exp(canon_vec)/(1-exp(canon_vec))^2)) <= 10))
})
