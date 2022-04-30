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

test_that("generate_data works for poisson", {
  trials <- 100
  set.seed(123)
  n <- 8
  p <- 8
  k <- 2
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- log(tcrossprod(x_mat, y_mat))
  library_size_vec <- 1:n

  ## Simulate data
  dat_array <- array(NA, dim = c(trials, n, p))
  for(i in 1:trials){
    set.seed(i)
    dat_array[i,,] <- generate_data(nat_mat, family = "poisson", library_size_vec = library_size_vec)
  }

  empirical_mean <- apply(dat_array, c(2,3), mean)
  empirical_var <- apply(dat_array, c(2,3), stats::var)
  empirical_mean2 <- apply(dat_array[1:100,,], c(2,3), mean)
  empirical_var2 <- apply(dat_array[1:100,,], c(2,3), stats::var)

  predicted_mean <- matrix(0, nrow = n, ncol = p)
  predicted_var <- matrix(0, nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      predicted_mean[i,j] <- library_size_vec[i]*exp(nat_mat[i,j])
      predicted_var[i,j] <- library_size_vec[i]*exp(nat_mat[i,j])
    }
  }

  expect_true(sum(abs(empirical_mean - predicted_mean)) <= sum(abs(empirical_mean2 - predicted_mean)))
  expect_true(sum(abs(empirical_var - predicted_var)) <= sum(abs(empirical_var2 - predicted_var)))

  # try many permutations
  error_vec <- c(sum(abs(empirical_mean - predicted_mean)/predicted_mean),
                 sum(abs(empirical_var - predicted_var)/predicted_var))
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    idx1 <- sample(1:n); idx2 <- sample(1:p)
    error_vec2 <- c(sum(abs(empirical_mean[idx1,idx2] - predicted_mean)/predicted_mean),
                    sum(abs(empirical_var[idx1,idx2] - predicted_var)/predicted_var))

    error_vec[1] <= error_vec2[1] | error_vec[2] <= error_vec2[2]
  })

  expect_true(all(bool_vec))
})

test_that("generate_data works for negative binomial", {
  trials <- 100
  set.seed(123)
  n <- 8
  p <- 8
  k <- 2
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  library_size_vec <- 1:n
  nuisance_param_vec <- c(1:p)*10

  ## Simulate data
  dat_array <- array(NA, dim = c(trials, n, p))
  for(i in 1:trials){
    set.seed(i)
    dat_array[i,,] <- eSVD2::generate_data(
      nat_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec, library_size_vec = library_size_vec
    )
  }

  empirical_mean <- apply(dat_array, c(2,3), mean)
  empirical_var <- apply(dat_array, c(2,3), stats::var)
  empirical_mean2 <- apply(dat_array[1:100,,], c(2,3), mean)
  empirical_var2 <- apply(dat_array[1:100,,], c(2,3), stats::var)


  predicted_mean <- matrix(0, nrow = n, ncol = p)
  predicted_var <- matrix(0, nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      predicted_mean[i,j] <- library_size_vec[i] * nuisance_param_vec[j] * exp(nat_mat[i,j])/(1-exp(nat_mat[i,j]))
      predicted_var[i,j] <- library_size_vec[i] * nuisance_param_vec[j] * exp(nat_mat[i,j])/(1-exp(nat_mat[i,j]))^2
    }
  }

  expect_true(sum(abs(empirical_mean - predicted_mean)) <= sum(abs(empirical_mean2 - predicted_mean)))
  expect_true(sum(abs(empirical_var - predicted_var)) <= sum(abs(empirical_var2 - predicted_var)))

  # try many permutations
  error_vec <- c(sum(abs(empirical_mean - predicted_mean)/predicted_mean),
                 sum(abs(empirical_var - predicted_var)/predicted_var))
  bool_vec <- sapply(1:trials, function(x){
    set.seed(x)
    idx1 <- sample(1:n); idx2 <- sample(1:p)
    error_vec2 <- c(sum(abs(empirical_mean[idx1,idx2] - predicted_mean)/predicted_mean),
                    sum(abs(empirical_var[idx1,idx2] - predicted_var)/predicted_var))

    error_vec[1] <= error_vec2[1] | error_vec[2] <= error_vec2[2]
  })

  expect_true(all(bool_vec))
})
