context("Test distribution")

## compute_mean is corrrect

test_that("compute_mean works", {
  mat <- matrix(1:25, 5, 5)
  res <- compute_mean(mat, family = "gaussian")

  expect_true(all(mat == res))
})

test_that("compute_mean is correct for gaussian", {
  trials <- 1000
  nat_mat <- matrix(c(1,10,20,30,40,100), nrow = 2, ncol = 3)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "gaussian", nuisance_param_vec = 1)
  }

  mean_mat <- apply(dat_array, c(1,2), mean)
  expected_mean_mat <- compute_mean(nat_mat, family = "gaussian")

  expect_true(all(abs(mean_mat - expected_mean_mat)/expected_mean_mat <= 0.1))
})

test_that("compute_mean is correct for curved gaussian", {
  trials <- 1000
  nat_mat <- matrix(1/c(1,10,20,30,40,100), nrow = 2, ncol = 3)
  nuisance_param_vec <- c(1,2,4)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "curved_gaussian",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(dat_array, c(1,2), mean)
  expected_mean_mat <- compute_mean(nat_mat, family = "curved_gaussian",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(mean_mat - expected_mean_mat)/expected_mean_mat <= 0.1))
})

test_that("compute_mean is correct for poisson", {
  trials <- 1000
  nat_mat <- matrix(log(c(1,10,20,30,40,100)), nrow = 2, ncol = 3)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "poisson",
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(dat_array, c(1,2), mean)
  expected_mean_mat <- compute_mean(nat_mat, family = "poisson",
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(mean_mat - expected_mean_mat)/expected_mean_mat <= 0.1))
})

test_that("compute_mean is correct for exponential", {
  trials <- 1000
  nat_mat <- matrix(-1/c(1,10,20,30,40,100), nrow = 2, ncol = 3)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "exponential",
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(dat_array, c(1,2), mean)
  expected_mean_mat <- compute_mean(nat_mat, family = "exponential",
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(mean_mat - expected_mean_mat)/expected_mean_mat <= 0.1))
})

test_that("compute_mean is correct for negative binomial", {
  trials <- 1000
  nat_mat <- matrix(log(c(0.9, 0.7, 0.5, 0.3, 0.1, 0.001)), nrow = 2, ncol = 3)
  nuisance_param_vec <- c(10, 50, 100)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "neg_binom",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)
  }

  mean_mat <- apply(dat_array, c(1,2), mean)
  expected_mean_mat <- compute_mean(nat_mat, family = "neg_binom",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(mean_mat - expected_mean_mat)/expected_mean_mat <= 0.1))
})

###########################################

## .compute_variance is correct


test_that(".compute_variance works", {
  mat <- matrix(log(1:25), 5, 5)
  res <- compute_mean(mat, family = "poisson")

  expect_true(all(dim(res) == dim(mat)))
})

test_that(".compute_variance is correct for gaussian", {
  trials <- 1000
  nat_mat <- matrix(c(1,10,20,30,40,100), nrow = 2, ncol = 3)
  nuisance_param_vec <- c(1,2,4)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "gaussian",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)
  }

  var_mat <- apply(dat_array, c(1,2), stats::var)
  expected_var_mat <- .compute_variance(nat_mat, family = "gaussian",
                                        nuisance_param_vec = nuisance_param_vec,
                                        library_size_vec = library_size_vec)

  expect_true(all(abs(var_mat - expected_var_mat)/expected_var_mat <= 0.1))
})

test_that(".compute_variance is correct for curved gaussian", {
  trials <- 1000
  nat_mat <- matrix(1/c(10,10,20,30,40,100), nrow = 2, ncol = 3)
  nuisance_param_vec <- c(2,4,10)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "curved_gaussian",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)
  }

  var_mat <- apply(dat_array, c(1,2), stats::var)
  expected_var_mat <- .compute_variance(nat_mat, family = "curved_gaussian",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(var_mat - expected_var_mat)/expected_var_mat <= 0.1))
})

test_that(".compute_variance is correct for poisson", {
  trials <- 1000
  nat_mat <- matrix(log(c(1,10,20,30,40,100)), nrow = 2, ncol = 3)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "poisson",
                                    library_size_vec = library_size_vec)
  }

  var_mat <- apply(dat_array, c(1,2), stats::var)
  expected_var_mat <- .compute_variance(nat_mat, family = "poisson",
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(var_mat - expected_var_mat)/expected_var_mat <= 0.1))
})

test_that(".compute_variance is correct for exponential", {
  trials <- 1000
  nat_mat <- matrix(-1/c(1,10,20,30,40,100), nrow = 2, ncol = 3)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "exponential",
                                    library_size_vec = library_size_vec)
  }

  var_mat <- apply(dat_array, c(1,2), stats::var)
  expected_var_mat <- .compute_variance(nat_mat, family = "exponential",
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(var_mat - expected_var_mat)/expected_var_mat <= 0.1))
})

test_that(".compute_variance is correct for negative binomial", {
  trials <- 1000
  nat_mat <- matrix(log(c(0.9, 0.7, 0.5, 0.3, 0.1, 0.001)), nrow = 2, ncol = 3)
  nuisance_param_vec <-  c(10, 50, 100)
  library_size_vec <- c(1,5)

  dat_array <- array(NA, c(2,3,trials))
  for(i in 1:trials){
    set.seed(i)
    dat_array[,,i] <- generate_data(nat_mat, family = "neg_binom",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)
  }

  var_mat <- apply(dat_array, c(1,2), stats::var)
  expected_var_mat <- .compute_variance(nat_mat, family = "neg_binom",
                                    nuisance_param_vec = nuisance_param_vec,
                                    library_size_vec = library_size_vec)

  expect_true(all(abs(var_mat - expected_var_mat)/expected_var_mat <= 0.1))
})
