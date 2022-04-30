context("Test initialization")

## initialize_esvd is correct

test_that("initialize_esvd works", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
  covariates <- matrix(rnorm(n * 3, mean = -1), nrow = n, ncol = 3)
  colnames(covariates) <- paste0("covariate_", 1:3)
  z_mat <- cbind(rep(1,p), rep(0.1,p), rep(.5,p))
  colnames(z_mat) <- paste0("covariate_", 1:3)
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)
  nat_mat <- pmin(nat_mat, 2)

  # Simulate data
  dat <- generate_data(exp(nat_mat),
                       family = "poisson",
                       nuisance_param_vec = NA,
                       library_size_vec = 1)

  res <- initialize_esvd(dat,
                         covariates = covariates,
                         case_control_variable = "covariate_1",
                         k = 2,
                         lambda = 0.1,
                         mixed_effect_variables = c("covariate_2", "covariate_3"))


  expect_true(is.list(res))
  expect_true(all(sort(names(res)) == sort(c("x_mat", "y_mat",
                                             "covariates", "z_mat",
                                             "log_pval_vec","param"))))
})
