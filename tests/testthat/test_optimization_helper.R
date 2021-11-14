context("Test optimization helper functions")

## .opt_nuisance is correct

test_that(".opt_nuisance works", {
  # Simulate data
  set.seed(123)
  n <- 300
  p <- 100
  k <- 2
  nuisance_param_vec <- runif(p, 0, 10)
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)

  dat <- eSVD2::generate_data(
    nat_mat, family = "neg_binom2", nuisance_param_vec = nuisance_param_vec,
    library_size_vec = 1
  )

  res <- .opt_nuisance(covariates = NULL,
                       dat = dat,
                       gene_group_factor = factor(rep("1", p)),
                       max_cell_subsample = 10*n,
                       offset_vec = rep(0, n),
                       x_mat = x_mat,
                       yb_mat = y_mat,
                       value_lower = NA,
                       value_upper = NA,
                       verbose = 0)
  expect_true(length(res) == p)
  expect_true(is.numeric(res))
  expect_true(all(res > 0))

  res <- .opt_nuisance(covariates = NULL,
                       dat = dat,
                       gene_group_factor = factor(1:p),
                       max_cell_subsample = 10*n,
                       offset_vec = rep(0, n),
                       x_mat = x_mat,
                       yb_mat = y_mat,
                       value_lower = rep(0.1, p),
                       value_upper = rep(NA, p),
                       verbose = 0)
  expect_true(length(res) == p)
  expect_true(is.numeric(res))
  expect_true(all(res > 0))
})

test_that(".opt_nuisance works for a more challenging setting", {
  set.seed(10)
  n <- 2000
  p <- 10
  k <- 2
  nuisance_param_vec <- runif(p, min = 0, max = 1)
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  nuisance_param_vec <- c(rep(0.1, 3), rep(1, 3), rep(500, 4))

  dat <- matrix(NA, n, p)
  for(i in 1:n){
    for(j in 1:p){
      dat[i,j] <- stats::rnbinom(n = 1, size = nuisance_param_vec[j],
                                 mu = exp(nat_mat[i,j]))
    }
  }

  res <- .opt_nuisance(covariates = NULL,
                       dat = dat,
                       gene_group_factor = factor(1:p),
                       max_cell_subsample = 10*n,
                       offset_vec = rep(0, n),
                       x_mat = x_mat,
                       yb_mat = y_mat,
                       value_lower = rep(0.1, p),
                       value_upper = rep(500, p),
                       verbose = 0)

  expect_true(sum(abs(res[1:3] - 0.1)) <= 1e-3)
  expect_true(sum(abs(res[4:6] - 1)) <= min(abs(res[-c(4:6)] - 1)))
  expect_true(sum(abs(res[7:10] - 500)) <= 1e-3)
})


