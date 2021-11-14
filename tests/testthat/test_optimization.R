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
                                             "offset_vec", "covariates", "param"))))
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
                                             "offset_vec", "covariates", "param"))))
})


test_that("opt_esvd works for negative binomial2, with different gene groups", {
  set.seed(10)
  n <- 500
  p <- 90
  k <- 2
  x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
  y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
  nat_mat <- tcrossprod(x_mat, y_mat)
  nuisance_param_vec <- c(runif(p/3, min = 0.05, max = 0.2),
                          runif(p/3, min = 0.9, max = 1.1),
                          runif(p/3, min = 300, max = 700))

  ## Simulate data
  dat <- matrix(NA, n, p)
  for(i in 1:n){
    for(j in 1:p){
      dat[i,j] <- stats::rnbinom(n = 1, size = nuisance_param_vec[j],
                                 mu = exp(nat_mat[i,j]))
    }
  }

  tmp <- log1p(dat)
  svd_res <- irlba::irlba(tmp, nv = 2)
  init_x_mat <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  init_y_mat <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  ## Optimization
  res <- opt_esvd(init_x_mat, init_y_mat, dat, family = "neg_binom2",
                  nuisance_param_vec = c(rep(0.1, p/3), rep(0.1, p/3), rep(500, p/3)),
                  library_size_vec = NA,
                  max_iter = 15,
                  gene_group_factor = factor(c(rep("1", p/3), rep("2", p/3), rep("3", p/3))),
                  nuisance_value_lower = c(0.05, 0.05, 1),
                  nuisance_value_upper = c(1, 1, 1000),
                  reestimate_nuisance = T,
                  reestimate_nuisance_per_iteration = 3,
                  verbose = 0)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat", "b_mat", "loss",
                                             "nuisance_param_vec", "library_size_vec",
                                             "offset_vec", "covariates", "param"))))
})
