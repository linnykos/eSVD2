context("Test gamma rate")

## gamma_rate is correct

test_that("gamma_rate and log_gamma_rate work", {
  set.seed(123)
  n <- 100
  p <- 150
  k <- 5
  x_mat <- matrix(abs(rnorm(n * k))*.5, nrow = n, ncol = k)
  y_mat <- matrix(abs(rnorm(p * k))*.5, nrow = p, ncol = k)
  covariates <- cbind(c(rep(0, n/2), rep(1, n/2)),
                      matrix(abs(rnorm(n * 3, mean = 1, sd = 0.1)), nrow = n, ncol = 3))
  colnames(covariates) <- paste0("covariate_", 1:4)
  z_mat <- cbind(c(rep(0, p/2), rep(2, p/2)), rep(1,p), rep(1,p), rep(1,p))
  colnames(z_mat) <-  colnames(covariates)
  case_control_variable <- "covariate_1"
  case_control_idx <- which(colnames(z_mat) == case_control_variable)
  nat_mat_nolib <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx])
  library_mat <- exp(tcrossprod(covariates[,-case_control_idx], z_mat[,-case_control_idx]))
  nuisance_vec <- rep(c(5, 1, 1/5), times = 50)

  # Simulate data
  gamma_mat <- matrix(NA, nrow = n, ncol = p)
  dat <- matrix(NA, nrow = n, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      gamma_mat[i,j] <- stats::rgamma(n = 1,
                                      shape = nuisance_vec[j]*exp(nat_mat_nolib[i,j]),
                                      rate = nuisance_vec[j])
      dat[i,j] <- stats::rpois(n = 1, lambda = library_mat[i,j] * gamma_mat[i,j])
    }
  }
  dat <- pmin(dat, 200)
  dat <- Matrix::Matrix(dat, sparse = T)
  rownames(dat) <- paste0("c", 1:n)
  colnames(dat) <- paste0("g", 1:p)
  metadata <- data.frame(individual = factor(rep(1:4, each = n/4)))
  rownames(metadata) <- rownames(dat)

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
                       max_iter = 10)

  nat_mat1 <- tcrossprod(eSVD_obj$fit_First$x_mat, eSVD_obj$fit_First$y_mat)
  nat_mat2 <- tcrossprod(eSVD_obj$covariates[,case_control_variable,drop = F],
                         eSVD_obj$fit_First$z_mat[,case_control_variable,drop = F])
  nat_mat_nolib <- nat_mat1 + nat_mat2
  mean_mat_nolib <- exp(nat_mat_nolib)

  offset_variable <- setdiff(colnames(eSVD_obj$covariates), case_control_variable)
  library_mat <- exp(tcrossprod(
    eSVD_obj$covariates[,offset_variable], eSVD_obj$fit_First$z_mat[,offset_variable]
  ))

  bool_vec <- sapply(1:ncol(dat), function(j){
    res1 <- gamma_rate(x = dat[,j],
                       mu = mean_mat_nolib[,j],
                       s = library_mat[,j])
    bool1 <- length(res1) == 1 & res1 >= 0

    res2 <- log_gamma_rate(x = dat[,j],
                           mu = mean_mat_nolib[,j],
                           s = library_mat[,j])
    bool2 <- length(res2) == 1 & exp(res2) >= 0

    bool3 <- abs(res1 - exp(res2)) <= 1e-3

    bool1 & bool2 & bool3
  })

  expect_true(all(bool_vec))
})
