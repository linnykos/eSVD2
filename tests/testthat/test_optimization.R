context("Test optimization")

test_that("opt_esvd works for eSVD_obj", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))/10, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))/10, nrow = p, ncol = k)
  covariates <- matrix(rnorm(n * 4, mean = 1, sd = 0.1), nrow = n, ncol = 4)
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(1, p/2)), rep(1,p), rep(1,p), rep(0.2,p))
  colnames(z_mat) <- paste0("covariate_", 1:4)
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)

  # Simulate data
  dat <- generate_data(nat_mat,
                       family = "poisson",
                       nuisance_param_vec = NA,
                       library_size_vec = 1)
  dat <- Matrix::Matrix(dat, sparse = T)

  eSVD_obj <- initialize_esvd(dat = dat,
                              covariates = covariates,
                              case_control_variable = "covariate_1",
                              k = 2,
                              lambda = 0.1,
                              mixed_effect_variables = c("covariate_2", "covariate_3"),
                              offset_variables = "covariate_3")
  eSVD_obj <- apply_initial_threshold(eSVD_obj = eSVD_obj,
                                      pval_thres = 0.1)
  res <- opt_esvd(input_obj = eSVD_obj,
                  max_iter = 5)

  expect_true(is.list(res))
  expect_true(inherits(res, "eSVD"))
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates",
                                             "initial_Reg", "param", "fit_Init",
                                             "fit_First", "latest_Fit"))))
  expect_true(all(sort(names(res$fit_First)) == sort(c("x_mat", "y_mat",
                                                      "z_mat", "loss"))))
  expect_true(all(diff(res$fit_First$loss) < 0))
  expect_true(res$latest_Fit == "fit_First")
})

test_that("opt_esvd works for sparse matrices", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))/10, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))/10, nrow = p, ncol = k)
  covariates <- matrix(rnorm(n * 4, mean = 1, sd = 0.1), nrow = n, ncol = 4)
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(1, p/2)), rep(1,p), rep(1,p), rep(0.2,p))
  colnames(z_mat) <- paste0("covariate_", 1:4)
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)

  # Simulate data
  dat <- generate_data(nat_mat,
                       family = "poisson",
                       nuisance_param_vec = NA,
                       library_size_vec = 1)
  dat <- Matrix::Matrix(dat, sparse = T)

  eSVD_obj <- initialize_esvd(dat = dat,
                              covariates = covariates,
                              case_control_variable = "covariate_1",
                              k = 2,
                              lambda = 0.1,
                              mixed_effect_variables = c("covariate_2", "covariate_3"),
                              offset_variables = "covariate_3")
  eSVD_obj <- apply_initial_threshold(eSVD_obj = eSVD_obj,
                                      pval_thres = 0.1)
  res <- opt_esvd(input_obj = eSVD_obj$dat,
                  x_init = eSVD_obj$fit_Init$x_mat,
                  y_init = eSVD_obj$fit_Init$y_mat,
                  z_init = eSVD_obj$fit_Init$z_mat,
                  covariates = eSVD_obj$covariates,
                  max_iter = 5)

  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat", "loss",
                                             "z_mat", "param", "covariates"))))
  expect_true(all(diff(res$loss) < 0))
})
