context("Test posterior")

## compute_posterior is correct

test_that("compute_posterior works", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))*.5, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))*.5, nrow = p, ncol = k)
  covariates <- cbind(c(rep(0, n/2), rep(1, n/2)),
                      matrix(rnorm(n * 3, mean = 1, sd = 0.1), nrow = n, ncol = 3))
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(1, p/2)), rep(1,p), rep(1,p), rep(0.2,p))
  colnames(z_mat) <-  colnames(covariates)
  case_control_variable <- "covariate_1"
  case_control_idx <- which(colnames(z_mat) == case_control_variable)
  nat_mat_nolib <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx])
  library_mat <- tcrossprod(covariates[,-case_control_idx], z_mat[,-case_control_idx])
  nuisance_vec <- rep(c(5, 1, 1/5), times = 50)

  # Simulate data
  gamma_mat <- matrix(NA, nrow = n, ncol = p)
  dat <- matrix(NA, nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      gamma_mat[i,j] <- stats::rgamma(n = 1,
                                      shape = nuisance_vec[j]*nat_mat_nolib[i,j],
                                      rate = nuisance_vec[j])
      dat[i,j] <- stats::rpois(n = 1, lambda = library_mat[i,j] * gamma_mat[i,j])
    }
  }
  dat <- Matrix::Matrix(dat, sparse = T)
  rownames(dat) <- paste0("c", 1:n)
  colnames(dat) <- paste0("g", 1:p)

  # fit eSVD
  eSVD_obj <- initialize_esvd(dat = dat,
                              covariates = covariates,
                              case_control_variable = case_control_variable,
                              k = 5,
                              lambda = 0.1,
                              mixed_effect_variables = c("covariate_2", "covariate_3", "covariate_4"),
                              offset_variables = NULL)
  eSVD_obj <- apply_initial_threshold(eSVD_obj = eSVD_obj,
                                      pval_thres = 0.1)
  eSVD_obj <- opt_esvd(input_obj = eSVD_obj,
                       max_iter = 5)
  eSVD_obj <- estimate_nuisance(input_obj = eSVD_obj,
                                verbose = 0)
  res <- compute_posterior(input_obj = eSVD_obj)

  expect_true(all(c("posterior_mean_mat", "posterior_var_mat") %in% names(res$fit_First)))
  expect_true(all(dim(res$fit_First$posterior_mean_mat) == dim(dat)))
  expect_true(all(dim(res$fit_First$posterior_var_mat) == dim(dat)))
})
