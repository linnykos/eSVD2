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

test_that("generate_data for gaussian corrrectly respects library size and nuisance", {
  trials <- 500
  canon_vec <- c(10,20,30)
  nat_mat <- sapply(canon_vec, function(x){rep(x, 2)})
  library_size_vec <- c(1:2)
  nuisance_param_vec <- c(0.01,5,20)

  res_array <- array(NA, dim = c(nrow(nat_mat), ncol(nat_mat), trials))
  for(i in 1:trials){
    set.seed(i)
    res_array[,,i] <- generate_data(nat_mat, family = "gaussian",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(res_array, c(1,2), mean)
  sd_mat <- apply(res_array, c(1,2), stats::sd)

  expected_mean_mat <- diag(library_size_vec) %*% compute_mean(nat_mat, family = "gaussian")
  expected_sd_mat <- diag(sqrt(library_size_vec)) %*% matrix(rep(nuisance_param_vec, 2), nrow = 2, ncol = 3, byrow = T)

  expect_true(all(abs(mean_mat - expected_mean_mat) < 1))
  expect_true(all(abs(sd_mat - expected_sd_mat) < 2))
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

test_that("generate_data for poisson corrrectly respects library size", {
  trials <- 1000
  canon_vec <- c(10,20,30)
  nat_mat <- sapply(log(canon_vec), function(x){rep(x, 2)})
  library_size_vec <- c(1:2)

  res_array <- array(NA, dim = c(nrow(nat_mat), ncol(nat_mat), trials))
  for(i in 1:trials){
    set.seed(i)
    res_array[,,i] <- generate_data(nat_mat, family = "poisson",
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(res_array, c(1,2), mean)
  var_mat <- apply(res_array, c(1,2), stats::var)

  expected_mean_mat <- diag(library_size_vec) %*% compute_mean(nat_mat, family = "poisson")
  expected_var_mat <- expected_mean_mat

  expect_true(all(abs(mean_mat - expected_mean_mat) < 1))
  expect_true(all(abs(var_mat - expected_var_mat) < 5))
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

test_that("generate_data for curved gaussian respects tol", {
 set.seed(5)
 nat_mat <- matrix(1/c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2, tol = 1e-3)

 expect_true(all(res > 0))

 set.seed(5)
 res <- generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2, tol = NA)

 expect_true(!all(res > 0))
})

test_that("generate_data for curved gaussian corrrectly respects library size and nuisance", {
  trials <- 1000
  canon_vec <- c(10,20,30)
  nat_mat <- sapply(1/canon_vec, function(x){rep(x, 2)})
  nuisance_param_vec <- c(2,4,10)
  library_size_vec <- c(1:2)

  res_array <- array(NA, dim = c(nrow(nat_mat), ncol(nat_mat), trials))
  for(i in 1:trials){
    set.seed(i)
    res_array[,,i] <- generate_data(nat_mat, family = "curved_gaussian",
                                    library_size_vec = library_size_vec,
                                    nuisance_param_vec = nuisance_param_vec)
  }

  mean_mat <- apply(res_array, c(1,2), mean)
  sd_mat <- apply(res_array, c(1,2), stats::sd)

  expected_mean_mat <- diag(library_size_vec) %*% compute_mean(nat_mat, family = "curved_gaussian")
  expected_sd_mat <- diag(sqrt(library_size_vec)) %*% compute_mean(nat_mat, family = "curved_gaussian") / rep(nuisance_param_vec, each = 2)

  expect_true(all(abs(mean_mat - expected_mean_mat) < 1))
  expect_true(all(abs(sd_mat - expected_sd_mat) < 2))
})


############

test_that("generate_data works for exponential", {
 set.seed(5)
 nat_mat <- matrix(-1/c(1:30), 5, 6)
 res <- generate_data(nat_mat, family = "exponential")

 expect_true(is.matrix(res))
 expect_true(all(dim(res) == dim(nat_mat)))
})

test_that("generate_data has the correct mean and variance for exponential", {
 set.seed(5)
 canon_vec <- seq(10, 50, length.out = 5)
 nat_mat <- sapply(-1/canon_vec, function(x){rep(x, 2000)})
 res <- generate_data(nat_mat, family = "exponential")

 mean_vec <- colMeans(res)
 expect_true(all(abs(mean_vec - canon_vec) <= 2))

 sd_vec <- apply(res, 2, stats::sd)
 expect_true(all(abs(sd_vec - canon_vec) <= 3))
})

test_that("generate_data for exponential corrrectly respects library size", {
  trials <- 1000
  canon_vec <- c(10,20,30)
  nat_mat <- sapply(-1/canon_vec, function(x){rep(x, 2)})
  library_size_vec <- c(1:2)

  res_array <- array(NA, dim = c(nrow(nat_mat), ncol(nat_mat), trials))
  for(i in 1:trials){
    set.seed(i)
    res_array[,,i] <- generate_data(nat_mat, family = "exponential",
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(res_array, c(1,2), mean)
  sd_mat <- apply(res_array, c(1,2), stats::sd)

  expected_mean_mat <- diag(library_size_vec) %*% compute_mean(nat_mat, family = "exponential")
  expected_sd_mat <- expected_mean_mat

  expect_true(all(abs(mean_mat - expected_mean_mat) < 2))
  expect_true(all(abs(sd_mat - expected_sd_mat) < 2))
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

test_that("generate_data for neg_binom corrrectly respects library size and nuisance", {
  trials <- 1000
  canon_vec <- c(10,20,30)
  nat_mat <- sapply(-1/canon_vec, function(x){rep(x, 2)})
  nuisance_param_vec <- c(10, 50, 100)
  library_size_vec <- c(1:2)

  res_array <- array(NA, dim = c(nrow(nat_mat), ncol(nat_mat), trials))
  for(i in 1:trials){
    set.seed(i)
    res_array[,,i] <- generate_data(nat_mat, family = "neg_binom",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(res_array, c(1,2), mean)
  var_mat <- apply(res_array, c(1,2), stats::var)

  expected_mean_mat <- diag(library_size_vec) %*% compute_mean(nat_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec)
  expected_var_mat <- expected_mean_mat/(1-.convert_natural_to_canonical(nat_mat, family = "neg_binom"))

  expect_true(all(abs(mean_mat - expected_mean_mat)/expected_mean_mat < 0.1))
  expect_true(all(abs(var_mat - expected_var_mat)/expected_var_mat < 0.1))
})
