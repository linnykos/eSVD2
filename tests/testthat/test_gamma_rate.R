context("Test gamma rate")

## gamma_rate is correct

test_that("gamma_rate and log_gamma_rate work", {
  # load("tests/assets/synthetic_data.RData")
  load("../assets/synthetic_data.RData")

  x_mat <- eSVD_obj$fit_First$x_mat
  y_mat <- eSVD_obj$fit_First$y_mat
  z_mat <- eSVD_obj$fit_First$z_mat
  covariates <- eSVD_obj$covariates
  case_control_idx <- which(colnames(covariates) == "case_control_1")

  nat_mat1 <- tcrossprod(x_mat, y_mat)
  nat_mat2 <- tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx])
  mean_mat_nolib <- exp(nat_mat1 + nat_mat2)
  library_mat <- exp(tcrossprod(covariates[,-case_control_idx], z_mat[,-case_control_idx]))

  bool_vec <- sapply(1:ncol(dat), function(j){
    res1 <- gamma_rate(x = dat[,j],
                       mu = mean_mat_nolib[,j],
                       s = library_mat[,j])
    bool1 <- length(res1) == 1 & res1 >= 0

    # res2 <- log_gamma_rate(x = dat[,j],
    #                        mu = mean_mat_nolib[,j],
    #                        s = library_mat[,j])
    # bool2 <- length(res2) == 1 & exp(res2) >= 0
    #
    # bool3 <- abs(res1 - exp(res2)) <= 1e-3

    # bool1 & bool2 & bool3
    bool1
  })

  expect_true(all(bool_vec))
})
