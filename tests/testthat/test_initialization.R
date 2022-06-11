context("Test initialization")

## initialize_esvd is correct

test_that("initialize_esvd works", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))/10, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))/10, nrow = p, ncol = k)
  covariates <- cbind(sample(c(0,1), size = n, replace = T),
                      matrix(rnorm(n * 3, mean = 1, sd = 0.1), nrow = n, ncol = 3))
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(1, p/2)), rep(1,p), rep(1,p), rep(0.2,p))
  colnames(z_mat) <- paste0("covariate_", 1:4)
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)

  # Simulate data
  dat <- generate_data(nat_mat,
                       family = "poisson",
                       nuisance_param_vec = NA,
                       library_size_vec = 1)
  colnames(dat) <- paste0("g", 1:ncol(dat))
  rownames(dat) <- paste0("c", 1:nrow(dat))

  res <- initialize_esvd(dat = dat,
                         covariates = covariates,
                         case_control_variable = "covariate_1",
                         k = 2,
                         lambda = 0.1,
                         mixed_effect_variables = c("covariate_2", "covariate_3"),
                         offset_variables = "covariate_3")

  expect_true(inherits(res, "eSVD"))
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates",
                                             "initial_Reg", "param"))))
  expect_true(sum(abs(res$dat - dat)) <= 1e-4)
  expect_true(sum(abs(res$covariates - cbind(1, covariates))) <= 1e-4)
  expect_true(all(sort(names(res$initial_Reg)) == sort(c("log_pval", "z_mat1", "z_mat2"))))
  expect_true(all(res$initial_Reg$log_pval <= 0))
  expect_true(length(names(res$initial_Reg$log_pval)) > 0)
  expect_true(all(names(res$initial_Reg$log_pval) == colnames(dat)))
  expect_true(mean(res$initial_Reg$log_pval[1:p/2]) >= mean(res$initial_Reg$log_pval[(p/2+1):p]))


  res <- initialize_esvd(dat = dat,
                         bool_intercept = F,
                         covariates = covariates,
                         case_control_variable = "covariate_1",
                         k = 2,
                         lambda = 0.1,
                         mixed_effect_variables = c("covariate_2", "covariate_3"),
                         offset_variables = "covariate_3")

  expect_true(inherits(res, "eSVD"))
  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates",
                                             "initial_Reg", "param"))))
})

###########################

## apply_initial_threshold is correct

test_that("apply_initial_threshold works", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))/10, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))/10, nrow = p, ncol = k)
  covariates <- cbind(sample(c(0,1), size = n, replace = T),
                      matrix(rnorm(n * 3, mean = 1, sd = 0.1), nrow = n, ncol = 3))
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(1, p/2)), rep(1,p), rep(1,p), rep(1,p))
  colnames(z_mat) <- paste0("covariate_", 1:4)
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)

  # Simulate data
  dat <- generate_data(nat_mat,
                       family = "poisson",
                       nuisance_param_vec = NA,
                       library_size_vec = 1)
  dat <- Matrix::Matrix(dat, sparse = T)
  colnames(dat) <- paste0("g", 1:ncol(dat))
  rownames(dat) <- paste0("c", 1:nrow(dat))

  eSVD_obj <- initialize_esvd(dat = dat,
                              covariates = covariates,
                              case_control_variable = "covariate_1",
                              k = 2,
                              lambda = 0.1,
                              mixed_effect_variables = c("covariate_2", "covariate_3"))
  res <- apply_initial_threshold(eSVD_obj = eSVD_obj,
                                 pval_thres = 0.1)

  expect_true(inherits(res, "eSVD"))
  expect_true(is.list(res))
  expect_true("init_pval_thres" %in% names(res$param))
  expect_true(res$param$init_pval_thres == 0.1)
  expect_true(all(sort(names(res)) == sort(c("dat", "covariates",
                                             "initial_Reg", "param", "fit_Init",
                                             "latest_Fit"))))
  expect_true(all(sort(names(res$fit_Init)) == sort(c("x_mat", "y_mat",
                                                      "z_mat"))))
  expect_true(all(dim(res$fit_Init$x_mat) == c(n,2)))
  expect_true(all(dim(res$fit_Init$y_mat) == c(p,2)))
  expect_true(all(dim(res$fit_Init$z_mat) == c(p,5)))
  expect_true(length(colnames(res$fit_Init$z_mat)) > 0)
  expect_true(all(colnames(res$fit_Init$z_mat) == c("Intercept", colnames(covariates))))
})

########################

## .initialize_coefficient is correct

test_that(".initialize_coefficient works for sparse matrices", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))/10, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))/10, nrow = p, ncol = k)
  covariates <- cbind(sample(c(0,1), size = n, replace = T),
                      matrix(rnorm(n * 3, mean = 1, sd = 0.1), nrow = n, ncol = 3))
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(100, p/2)), rep(10,p), rep(10,p), rep(1,p))
  colnames(z_mat) <- paste0("covariate_", 1:4)
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)/50

  # Simulate data
  dat <- generate_data(nat_mat,
                       family = "poisson",
                       nuisance_param_vec = NA,
                       library_size_vec = 1)
  dat <- Matrix::Matrix(dat, sparse = T)
  colnames(dat) <- paste0("g", 1:ncol(dat))
  rownames(dat) <- paste0("c", 1:nrow(dat))

  res <- .initialize_coefficient(bool_intercept = T,
                                 case_control_variable = "covariate_1",
                                 covariates = covariates,
                                 dat = dat,
                                 lambda = 0.1,
                                 mixed_effect_variables = c("covariate_2", "covariate_3"),
                                 offset_variables = NULL)

  expect_true(inherits(res, "initial_Reg"))
  expect_true(all(res$log_pval <= 0))
  expect_true(all(sort(names(res)) == sort(c("log_pval", "z_mat1", "z_mat2"))))
  expect_true(mean(res$log_pval[1:p/2]) >= mean(res$log_pval[(p/2+1):p]))
})
